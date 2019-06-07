
/*
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <iostream>
#include <seqan/align.h>
#include <seqan/bam_io.h>
#include<iostream>
#include<fstream>
#include <seqan/find.h>
#include <stdlib.h>
#include <time.h>
// #include <map>
*/

#include <seqan/seq_io.h>
#include <seqan/bam_io.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/align.h>
#include <seqan/find.h>
// #include <seqan/index.h>
#include <seqan/arg_parse.h>

using namespace seqan;
using namespace std;





int main(int argc, char const * argv[])
{
    ArgumentParser parser("Split Reads");

//     addOption(parser, ArgParseOption("b", "bam", "Path to the bam file", ArgParseArgument::INPUT_FILE, "IN"));
//     setRequired(parser, "bam");

    addOption(parser, ArgParseOption("f", "flexbarResult", "Path to the flexbar result (fasta)", ArgParseArgument::INPUT_FILE, "IN"));
    setRequired(parser, "flexbarResult");
//     setRequired(parser, "reads1");

    addOption(parser, ArgParseOption("o", "output", "Path to output files prefix", ArgParseArgument::OUTPUT_FILE, "OUT"));
//     setRequired(parser, "output");

    addOption(parser, ArgParseOption("l", "barcodeLength", "Shorting reads to this length", ArgParseArgument::INTEGER, "INT"));
//     addOption(parser, ArgParseOption("bL", "barcodeL", "Number of reads with are loaded at the same Time. Default = 100000", ArgParseArgument::INTEGER, "INT"));
//     setRequired(parser, "barcodeL");
//     addOption(parser, ArgParseOption("bS", "batchSize", "Number of reads with are loaded at the same Time. Default = 100000", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("v", "verbose", ""));

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    CharString dictPath, flexPath, bamPath, outputPath;
    int barcodeLength = 0;
//     int batchSize1 = 100000, barcodeLength;

//     getOptionValue(bamPath, parser, "bam");
    getOptionValue(flexPath, parser, "flexbarResult");
    getOptionValue(outputPath, parser, "output");
    getOptionValue(barcodeLength, parser, "barcodeLength");
//     getOptionValue(batchSize1, parser, "batchSize");
    bool verbose = isSet(parser, "verbose");



    //prepare flex result fasta
    StringSet<CharString> idsReads;
    SeqFileIn seqFileInFlex(toCString(flexPath));


    CharString outputRevcomp = outputPath;
    outputPath += "_left_tail_trimmed.fasta";
    outputRevcomp += "_right_tail_trimmed.fasta";

    SeqFileOut seqFileOut(toCString(outputPath));
    SeqFileOut seqFileOutRevcomp(toCString(outputRevcomp));

    while (!atEnd(seqFileInFlex))
    {
        Dna5String read;
        CharString id;
        try
        {

            readRecord(id, read, seqFileInFlex);

            Finder<CharString> finder(id);
            Pattern<CharString, Horspool> pattern("_Flexbar_removal_");
            find(finder, pattern);
            int test = beginPosition(finder);
            if(test > 0){
                Finder<CharString> finder(id);
                Pattern<CharString, Horspool> pattern(" revcomp");
                find(finder, pattern);
                int end = beginPosition(finder);
    //             std::cout << "end: " << end << "\n";i
                bool doRC = end > 0;
                if (doRC){
                    if (barcodeLength > 0 && length(read) > (barcodeLength + 1))
                        read = suffix(read, length(read) - barcodeLength - 1);
                    writeRecord(seqFileOutRevcomp, id, read);
                }else{
                    if (barcodeLength > 0 && length(read) > (barcodeLength + 1))
                        read = prefix(read, length(read) - barcodeLength - 1);
                    writeRecord(seqFileOut, id, read);
                }
            }
        }
        catch (IOError const & e)
        {
            std::cerr << "ERROR: could not copy record. " << e.what() << "\n";
        }
    }

    close(seqFileInFlex);


    std::cout << "Finished!\n";

    return 0;
}



