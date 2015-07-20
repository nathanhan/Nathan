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
  bw.retrieveIndex(dbindex);
  SnowTools::BamReadVector brv;

  while (contigparser.GetNextRead(contig,isreadvalid)) {
    cout << "outside " << contig.Qname() << " " << contig.Position() << " " << contig.Sequence() << endl;
    if (isreadvalid) {
      // BamReadVector brv;
      //int dum = 0;
      bw.alignSingleSequence(contig.Sequence(),contig.Qname(),contig,false);
      cout << "inside " << contig.Qname() << " " << contig.Position() << " "  << contig.Sequence() << endl;
     
    }

  }

  return 0;
}
