#include "SnowTools/SnowUtils.h"
#include "SnowTools/BamWalker.h"
#include "SnowTools/BamRead.h"
#include "SnowTools/BWAWrapper.h"
#include <iostream> 
#include <string>
#include <vector>
#include <ctime>
#include <unordered_map>
#include <fstream>

// TODO
// 2) look up std::unordered_map. Use this to keep track of number of alignments per virus, etc
// 3) find way to filter secondary alignments so you only keep good ones. Map Q? Length of match? NM tag? etc
// 5) command line options
// 6) build index bac
// 8) SW, BLAST integration
// 9) error management: if no cigar, no seq, etc

// misc stuff
// small test contig file: /home/unix/nhan/snow/finalfixed.tttt.contigs.sort.bam and finalbig.contigs.bam

using namespace std;

const time_t ctt = time(0);

int main(int argc, char** argv) {
  
  //get time and print it to verify latest version
  cout << "test combine " << asctime(localtime(&ctt)) << endl;

  //set up BamWalker input reader dependencies
  string contigBamPath = "/home/unix/jwala/tmp.contigs.bam"; //input file path here, contains contigs to search for
  SnowTools::BamWalker contigBamParser(contigBamPath); //construct BamWalker object to extract contigs from bam
  SnowTools::BamRead contig; //construct container to store extracted contigs
  bool isReadValid; //is the contig I just read in not messed up?
 
  //set up BWA alignment algorithm and BamWalker output writer dependencies
  //declare a few alignment function arguments ahead of time
  bool areSecondaryAlignsKept = false; //yes I should keep less likely alignments
  SnowTools::BamReadVector refseqSearchResultTempContainer, virusesSearchResultTempContainer, TEeukSearchResultTempContainer, TEallSearchResultTempContainer; //store alignments here
  //fetch refseq db
  SnowTools::BWAWrapper searchRefseq; //construct BWA object
  string refseqPath = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta"; //file path to db
  searchRefseq.retrieveIndex(refseqPath); //tell BWA object to load this path's index in anticipation of alignment job
  //init refseq results filewriter
  bam_hdr_t * refseq_header = searchRefseq.HeaderFromIndex(); //construct bam header from index and store in container object
  SnowTools::BamWalker refseq_writer; //construct output bam filewriter object
  refseq_writer.SetWriteHeader(refseq_header); //feed bam header into filewriter
  refseq_writer.OpenWriteBam("/xchip/gistic/Nathan/output_writer_files/realignmentToReference.bam"); //physically open output file stream
  //fetch viruses db
  SnowTools::BWAWrapper searchViruses;
  string virusesPath = "/home/unix/jwala/SnowmanFilter/viral.1.1.genomic.fna";
  searchViruses.retrieveIndex(virusesPath);
  //init viruses results filewriter
  bam_hdr_t * viruses_header = searchViruses.HeaderFromIndex();
  SnowTools::BamWalker viruses_writer;
  viruses_writer.SetWriteHeader(viruses_header);
  viruses_writer.OpenWriteBam("/xchip/gistic/Nathan/output_writer_files/realignmentToViruses.bam");
  //fetch repbase TE db (all euks)
  SnowTools::BWAWrapper searchTEeuk;
  string TEeukPath = "/xchip/gistic/Nathan/Repbase/repbase_euk.fasta";
  searchTEeuk.retrieveIndex(TEeukPath);
  //init repbase TE results filewriter
  bam_hdr_t * TEeuk_header = searchTEeuk.HeaderFromIndex();
  SnowTools::BamWalker TEeuk_writer;
  TEeuk_writer.SetWriteHeader(TEeuk_header);
  TEeuk_writer.OpenWriteBam("/xchip/gistic/Nathan/output_writer_files/realignmentToTEeuk.bam");
  //fetch repeatmasker TE db (humans only)
  SnowTools::BWAWrapper searchTEall;
  string TEallPath = "/home/unix/jwala/SnowmanFilter/repeatmasker/RepeatMasker/Libraries/20110419/homo_sapiens/all.rep.fa";
  searchTEall.retrieveIndex(TEallPath);
  //init repeatmasker TE results filewriter
  bam_hdr_t * TEall_header = searchTEall.HeaderFromIndex();
  SnowTools::BamWalker TEall_writer;
  TEall_writer.SetWriteHeader(TEall_header);
  TEall_writer.OpenWriteBam("/xchip/gistic/Nathan/output_writer_files/realignmentToTEall.bam");

  //declare CIGAR parsing and clipped sequence extraction dependencies ahead of time
  SnowTools::Cigar Cig; //for storing CIGAR string from contig
  string entireContigSeq, clippedContigSeq; //for storing complete contig sequence, for storing clipped segment only
  int positionCounter; //keeps track of CIGAR lengths cumulatively
  string origContigQname, modifiedQname; //for storing original name of contig, for storing auto-generated name of child segment
  int QnameCounter; //appendix used to generate name of child segment logically

  //open filewriter to store clipped segments before further processing in FASTA data file
  remove("/xchip/gistic/Nathan/output_writer_files/clips.fasta"); //kill file if it exists already
  ofstream clippedContigSeq_writer;
  clippedContigSeq_writer.open("/xchip/gistic/Nathan/output_writer_files/clips.fasta"); //create new version of now-deleted file

  //init progress updates
  std::cerr << "Starting contig parsing from bam..."  << std::endl;
  int readcounter = 0;
  int validreads = 0;
  //

  //actual processing starts here
  //perform the BAMWalk read-in of next contig
  while (contigBamParser.GetNextRead(contig,isReadValid)) {
    
    ++readcounter;

    //print progress update
    if (readcounter % 50 == 0)
      std::cerr << "...working on read " << readcounter << " which is " << contig << " of total of " << validreads << " valid reads so far " << std::endl;
    //

    if (isReadValid) {

      ++validreads;

      //store this contig's raw sequence, qname, and CIGAR into containers
      entireContigSeq = contig.Sequence();
      origContigQname = contig.Qname();
      Cig = contig.GetCigar();
      //reset some counters
      positionCounter = 0;
      QnameCounter = 1;

      //for every CIGAR phrase in CIGAR
      for (int i = 0; i < Cig.size(); i++) {

        /*need to deal with 25M1D25S here*/
        //if it's a softclip
        if (Cig[i].Type == 'S') {
          
          //extract clipped sequence segments from contig
          clippedContigSeq = entireContigSeq.substr(positionCounter, Cig[i].Length);
          positionCounter+=Cig[i].Length;
          modifiedQname = origContigQname.append('.' + to_string(QnameCounter));
          QnameCounter++;

          //add clipped sequence segments to intermediate FASTA data file
          clippedContigSeq_writer << ">" << modifiedQname << endl << clippedContigSeq << endl << endl;

          //look them up in DBs with BWA
          searchRefseq.alignSingleSequence(clippedContigSeq, modifiedQname, refseqSearchResultTempContainer, areSecondaryAlignsKept);
          searchViruses.alignSingleSequence(clippedContigSeq, modifiedQname, virusesSearchResultTempContainer, areSecondaryAlignsKept);
          searchTEeuk.alignSingleSequence(clippedContigSeq, modifiedQname, TEeukSearchResultTempContainer, areSecondaryAlignsKept);
          searchTEall.alignSingleSequence(clippedContigSeq, modifiedQname, TEallSearchResultTempContainer, areSecondaryAlignsKept);

          //write alignment results to bam files
          if (refseqSearchResultTempContainer.size() > 0) {
            for (auto& result : refseqSearchResultTempContainer) {
              result.AddIntTag("RA", refseqSearchResultTempContainer.size());
              refseq_writer.WriteAlignment(result); // pass this a BamRead
            }
          }
          if (refseqSearchResultTempContainer.size() > 0) {  
            for (auto& result : virusesSearchResultTempContainer) {
              result.AddIntTag("RA", virusesSearchResultTempContainer.size());
              viruses_writer.WriteAlignment(result); // pass this a BamRead
            }
          }
          if (refseqSearchResultTempContainer.size() > 0) {
            for (auto& result : TEeukSearchResultTempContainer) {
              result.AddIntTag("RA", TEeukSearchResultTempContainer.size());
              TEeuk_writer.WriteAlignment(result); // pass this a BamRead
            }
          }
          if (refseqSearchResultTempContainer.size() > 0) {
            for (auto& result : TEallSearchResultTempContainer) {
              result.AddIntTag("RA", TEallSearchResultTempContainer.size());
              TEall_writer.WriteAlignment(result); // pass this a BamRead
            }
          }

        }

        else if (Cig[i].Type != 'D') {

          positionCounter+=Cig[i].Length;

        }

      }

  // example of getting string of "chr" or virus etc that query was aligned to
  // check for segfault below
  // std::string chr_name = std::string(reference_header->target_name[ refseqSearchResultTempContainer[0].ChrID()]);
  
      }

    }

  clippedContigSeq_writer.close();
  return 0;

}
