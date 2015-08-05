#include "SnowTools/SnowUtils.h"
#include "SnowTools/BamWalker.h"
#include "SnowTools/BamRead.h"
#include "SnowTools/BWAWrapper.h"
#include <iostream> 
#include <string>
#include <vector>

// TODO
// 1) Write the alignments to BAM (viral, etc)
// 2) look up std::unordered_map. Use this to keep track of number of alignments per virus, etc
// 3) find way to filter secondary alignments so you only keep good ones. Map Q? Length of match? NM tag? etc

using namespace std;

int main(int argc, char** argv) {
  
  //verify no funky compile issues
  cout << "test Thur 2:28pm" << endl;

  //set up some dependencies for BamWalker raw data parsing
  //string contigBamPath = "/home/unix/nhan/fixed.tttt.contigs.sort.bam"; //input file path here, contains contigs to search for
  string contigBamPath = "/home/unix/jwala/tmp.contigs.bam"; //input file path here, contains contigs to search for
  SnowTools::BamWalker contigBamParser(contigBamPath); //create structure to extract contigs from bam
  SnowTools::BamRead contig; //create structure to store extracted contigs
  bool isReadValid; //ie is the contig I just read in not messed up?

  //set up some dependencies for BWA search/alignment
  SnowTools::BWAWrapper searchRefseq; //create structure to prep + execute search on refseq                                                       
  string refseqPath = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta"; //file path to refseq
  searchRefseq.retrieveIndex(refseqPath); //load-in refseq index to prepare for search
  
  SnowTools::BWAWrapper searchViruses; //create structure to prep + execute search on viruses
  string virusesPath = "/home/unix/jwala/SnowmanFilter/viral.1.1.genomic.fna"; //file path to virus db
  searchViruses.retrieveIndex(virusesPath); //load-in viruses index to prepare for search

  //SnowTools::BWAWrapper searchBac; //create structure to prep + execute search on bacteria
  //string bacPath = "/home/unix/jwala/SnowmanFilter/viral.bacteria.fa"; //file path to bacteria db
  //searchBac.retrieveIndex(bacPath); //load-in bacteria index to prepare for search
  
  //SnowTools::BWAWrapper searchTE; //transposable elements, part 1
  //string TEPath = "/home/unix/jwala/SnowmanFilter/repeatmasker/RepeatMasker/Libraries/20110419/homo_sapiens/longlib";
  //searchTE.retrieveIndex(TEPath);

  // make BAM headers for database
  bam_hdr_t * reference_header = searchRefseq.HeaderFromIndex();
  //bam_hdr_t * virus_header = searchViruses.HeaderFromIndex();
  
  // make a bam walker for writing new alignments
  SnowTools::BamWalker my_writer;
  my_writer.SetWriteHeader(reference_header);
  my_writer.OpenWriteBam("realignmentToReference.bam");
  
  bool areSecondaryAlignsKept = true; //ie yes I should keep worse matches of search

  //set up some dependencies for CIGAR-specific parsing
  //for storing clipped seqs:
  SnowTools::Cigar Cig;
  string entireContigSeq;
  vector<string> clippedContigSeqs;
  int positionCounter;
  //for storing Qnames in parallel:
  string origContigQname, modifiedQname;
  vector<string> Qnames;
  int QnameCounter;

  //set up some containers for actual BWA search
  SnowTools::BamReadClusterVector refseqSearchResultContainer, virusesSearchResultContainer, bacSearchResultContainer;
  SnowTools::BamReadVector refseqSearchResultTempContainer, virusesSearchResultTempContainer, bacSearchResultTempContainer;

  std::cerr << "...Staring contig parsing from bam"  << std::endl;
  int counter = 0;
  int read_valid = 0;

  //actual program starts here
  //perform the BAMWalk read-in of next contig
  while (contigBamParser.GetNextRead(contig,isReadValid)) {
    
    if (++counter % 10 == 0)
      std::cerr << "...working on read " << counter << " which is " << contig << " of total of " << read_valid << " valid reads so far " << std::endl;
    if (isReadValid)
      ++read_valid;
    

    if (isReadValid) {

	  //store this contig's sequence, qname, and CIGAR into variables
	  entireContigSeq = contig.Sequence();
	  origContigQname = contig.Qname();
      SnowTools::Cigar Cig = contig.GetCigar();
      //reset some counters
	  positionCounter = 0;//cumulative length tracker (needed since CIG lengths are split ex 129M29S)
	  QnameCounter = 0;//generates appendices to keep modified Qnames unique


      //extract clipped sequence segments from contig, Qnames in parallel
      //store both in string vectors
      for (int i = 0; i < Cig.size(); i++) {

      	if (Cig[i].Type == 'S') {
      		clippedContigSeqs.push_back(entireContigSeq.substr(positionCounter, Cig[i].Length));
      		positionCounter+=Cig[i].Length;

      		modifiedQname = origContigQname.append(to_string(QnameCounter));
      		Qnames.push_back(modifiedQname);
      		QnameCounter++;
      	}

      	else {
      		positionCounter+=Cig[i].Length;
      	}

      }
      
      //perform BWA search on databases
      //store results in second permanent container
      for (int i = 0; i < clippedContigSeqs.size(); i++) {

      	//search on refseq
      	searchRefseq.alignSingleSequence(clippedContigSeqs[i],Qnames[i],refseqSearchResultTempContainer, areSecondaryAlignsKept);
      	//refseqSearchResultContainer.push_back(refseqSearchResultTempContainer);
      	//search on viruses
      	searchViruses.alignSingleSequence(clippedContigSeqs[i],Qnames[i],virusesSearchResultTempContainer,areSecondaryAlignsKept);
      	//virusesSearchResultContainer.push_back(virusesSearchResultTempContainer);
      	//search on bacteria
      	//searchBac.alignSingleSequence(clippedContigSeqs[i],Qnames[i],bacSearchResultTempContainer,areSecondaryAlignsKept);
      	//bacSearchResultContainer.push_back(bacSearchResultTempContainer);

	// 
	

	// example of getting string of "chr" or virus etc that query was aligned to
	// check for segfault below
	// std::string chr_name = std::string(reference_header->target_name[ refseqSearchResultTempContainer[0].ChrID()]);
	
	// to write the alignments out to a bam file
	for (auto& i : refseqSearchResultTempContainer) {
	  i.AddIntTag("RA", refseqSearchResultTempContainer.size());
	  //i.AddIntTag("VA", refseqSearchResultTempContainer.size());
	  my_writer.WriteAlignment(i); // pass this a BamRead
	  break;
	}

      }

    }

  }

  return 0;

}
