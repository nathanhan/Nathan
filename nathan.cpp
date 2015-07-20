#include "SnowTools/SnowUtils.h"
#include "SnowTools/BamWalker.h"
#include "SnowTools/BamRead.h"
#include "SnowTools/BWAWrapper.h"
#include <iostream> 
#include <string>

using namespace std;

int main(int argc, char** argv) {
  cout << "test" << endl;
  string inputbampath = "/home/unix/nhan/fixed.tttt.contigs.sort.bam";
  SnowTools::BamRead contig;
  SnowTools::BamWalker contigparser(inputbampath);
  bool isreadvalid;
  SnowTools::BWAWrapper bw;
  //SnowTools::USeqVector v = {{getContigName(), getSequence()}};                                                       
  string dbindex = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta";
  string dbindex_viral = "/xchip/gistic/Jeremiah/Projects/SnowmanFilters/viral.1.1.genomic.fna";
  bw.retrieveIndex(dbindex);
  SnowTools::BamReadVector brv;
  SnowTools::BamRead x;

  while (contigparser.GetNextRead(contig,isreadvalid)) {
    cout << "outside " << contig.Qname() << " " << contig.Position() << " " << contig.Sequence() << endl;
    if (isreadvalid) {
      // BamReadVector brv;
      //int dum = 0;
      BamReadVector this_brv;
      bw.alignSingleSequence(contig.Sequence(),contig.Qname(),this_brv,false);
      brv.insert(brv.end(), this_brv.begin(), this_brv.end());
      x = brv.back();
      cout << "inside " << x.Qname() << " " << x.Position() << " "  << x.Sequence() << endl;
     
    }

  }

  return 0;
}
