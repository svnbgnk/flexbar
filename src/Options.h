/*
 *   Options.h
 *
 *   Author: jtr
 */

#ifndef FLEXBAR_OPTIONS_H
#define FLEXBAR_OPTIONS_H

#include <seqan/arg_parse.h>

#include "FlexbarIO.h"


struct Options{

	std::string readsFile, readsFile2, barReadsFile;
	std::string barcodeFile, adapterFile, barcode2File, adapter2File;
	std::string adapterSeq, targetName, logAlignStr, outCompression;

	bool isPaired, useAdapterFile, useNumberTag, useRemovalTag, randTag, logStdout;
	bool switch2Fasta, writeUnassigned, writeSingleReads, writeSingleReadsP, writeLengthDist;
	bool useStdin, useStdout, relaxRegion, revCompAdapter, qtrimPostRm, bNoMBV;
        bool logEverything;
        int prefix;

	int cutLen_begin, cutLen_end, cutLen_read, a_tail_len, b_tail_len;
	int qtrimThresh, qtrimWinSize, a_overhang;
	int maxUncalled, min_readLen, a_min_overlap, b_min_overlap, nThreads, bundleSize;
	int match, mismatch, gapCost, b_match, b_mismatch, b_gapCost;

	float a_errorRate, b_errorRate;

	flexbar::TrimEnd         end, b_end;
	flexbar::FileFormat      format;
	flexbar::QualityType     qual;
	flexbar::QualTrimType    qTrim;
	flexbar::LogAlign        logAlign;
	flexbar::CompressionType cmprsType;
	flexbar::RunType         runType;
	flexbar::BarcodeDetect   barDetect;
	flexbar::AdapterRemoval  adapRm;

	tbb::concurrent_vector<flexbar::TBar> barcodes, adapters, barcodes2, adapters2;

	std::ostream *out;
	std::fstream fstrmOut;

	Options(){
		readsFile      = "";
		readsFile2     = "";
		barReadsFile   = "";
		barcodeFile    = "";
		adapterFile    = "";
		barcode2File   = "";
		adapter2File   = "";
		outCompression = "";

		isPaired          = false;
		useAdapterFile    = false;
		useNumberTag      = false;
		useRemovalTag     = false;
		writeUnassigned   = false;
                logEverything     = false;
		writeSingleReads  = false;
		writeSingleReadsP = false;
		writeLengthDist   = false;
		switch2Fasta      = false;
		logStdout         = false;
		randTag           = false;
		useStdin          = false;
		useStdout         = false;
		relaxRegion       = false;
		revCompAdapter    = false;
		qtrimPostRm       = false;
		bNoMBV            = false;

                prefix        = 0;
		cutLen_begin  = 0;
		cutLen_end    = 0;
		cutLen_read   = 0;
		qtrimThresh   = 0;
		qtrimWinSize  = 0;
		a_tail_len    = 0;
		b_tail_len    = 0;
		b_min_overlap = 0;

		format    = flexbar::FASTA;
		qual      = flexbar::SANGER;
		qTrim     = flexbar::QOFF;
		logAlign  = flexbar::NONE;
		cmprsType = flexbar::UNCOMPRESSED;
		barDetect = flexbar::BOFF;
		adapRm    = flexbar::AOFF;
    }
};


const std::string getFlexbarBanner(const seqan::CharString version){

	std::string banner = "";

	banner += "               ________          __              \n";
	banner += "              / ____/ /__  _  __/ /_  ____ ______\n";
	banner += "             / /_  / / _ \\| |/ / __ \\/ __ `/ ___/\n";
	banner += "            / __/ / /  __/>  </ /_/ / /_/ / /    \n";
	banner += "           /_/   /_/\\___/_/|_/_.___/\\__._/_/     \n\n";

	banner += "Flexbar - flexible barcode and adapter removal, version ";

	append(banner, version);

	banner += "\nDeveloped with SeqAn, the library for sequence analysis\n";

	return banner;
}


const std::string getFlexbarCitation(){
	return "Matthias Dodt, Johannes T. Roehr, Rina Ahmed, Christoph Dieterich:\nFlexbar - flexible barcode and adapter processing for next-generation\nsequencing platforms. Biology 2012, 1(3):895-905.\n";
}


const std::string getFlexbarURL(){
	return "Available on github.com/seqan/flexbar\n";
}


void defineOptions(seqan::ArgumentParser &parser, const std::string version, const std::string date){

	using namespace seqan;

	typedef ArgParseArgument ARG;

	setVersion(parser, version);
	setDate(parser, date);

	// setCitation(parser, "\n\n" + getFlexbarCitation());

	// setAppName(parser, "");
	// setShortCopyright(parser, "");
	// setLongCopyright(parser, "");
	// ARG::OUTPUTPREFIX

	setShortDescription(parser, "flexible barcode and adapter removal");

	addUsageLine(parser, "\\fB-r\\fP reads [\\fB-b\\fP barcodes] [\\fB-a\\fP adapters] [options]");

	addOption(parser, ArgParseOption("hm", "man-help", "Print advanced options as man document."));
	addOption(parser, ArgParseOption("v", "versions", "Print Flexbar and SeqAn version numbers."));
	addOption(parser, ArgParseOption("c", "cite", "Show program reference for citation."));

	addSection(parser, "Basic options");
	addOption(parser, ArgParseOption("n", "threads", "Number of threads to employ.", ARG::INTEGER));
	addOption(parser, ArgParseOption("N", "bundle", "Number of read pairs per thread.", ARG::INTEGER));
	addOption(parser, ArgParseOption("t", "target", "Prefix for output file names or paths.", ARG::STRING));
	addOption(parser, ArgParseOption("r", "reads", "Fasta/q file or stdin (-) with reads that may contain barcodes.", ARG::INPUT_FILE));
	addOption(parser, ArgParseOption("p", "reads2", "Second input file of paired reads, gz and bz2 files supported.", ARG::INPUT_FILE));

	addSection(parser, "Barcode detection");
	addOption(parser, ArgParseOption("b",  "barcodes", "Fasta file with barcodes for demultiplexing, may contain N.", ARG::INPUT_FILE));
	addOption(parser, ArgParseOption("b2", "barcodes2", "Additional barcodes file for second read set in paired mode.", ARG::INPUT_FILE));
	addOption(parser, ArgParseOption("br", "barcode-reads", "Fasta/q file containing separate barcode reads for detection.", ARG::INPUT_FILE));
	addOption(parser, ArgParseOption("be", "barcode-trim-end", "Type of detection, see section trim-end modes.", ARG::STRING));
	addOption(parser, ArgParseOption("bn", "barcode-tail-length", "Region size in tail trim-end modes. Default: barcode length.", ARG::INTEGER));
	addOption(parser, ArgParseOption("bo", "barcode-min-overlap", "Minimum overlap of barcode and read. Default: barcode length.", ARG::INTEGER));
	addOption(parser, ArgParseOption("bt", "barcode-error-rate", "Error rate threshold for mismatches and gaps.", ARG::DOUBLE));
	addOption(parser, ArgParseOption("bk", "barcode-keep", "Keep barcodes within reads instead of removal."));
	addOption(parser, ArgParseOption("bu", "barcode-unassigned", "Include unassigned reads in output generation."));
	// addOption(parser, ArgParseOption("ba", "barcode-no-mbv", "Turn off bit-vector alignment with edit distance."));
	addOption(parser, ArgParseOption("bm", "barcode-match", "Alignment match score.", ARG::INTEGER));
	addOption(parser, ArgParseOption("bi", "barcode-mismatch", "Alignment mismatch score.", ARG::INTEGER));
	addOption(parser, ArgParseOption("bg", "barcode-gap", "Alignment gap score.", ARG::INTEGER));
	// addOption(parser, ArgParseOption("bm", "barcode-match", "Alignment match score in case of no-mbv.", ARG::INTEGER));
	// addOption(parser, ArgParseOption("bi", "barcode-mismatch", "Standard alignment mismatch score.", ARG::INTEGER));
	// addOption(parser, ArgParseOption("bg", "barcode-gap", "Standard alignment gap score.", ARG::INTEGER));

	addSection(parser, "Adapter removal");
	addOption(parser, ArgParseOption("a",  "adapters", "Fasta file with adapters for removal that may contain N.", ARG::INPUT_FILE));
	addOption(parser, ArgParseOption("a2", "adapters2", "File with extra adapters for second read set in paired mode.", ARG::INPUT_FILE));
	addOption(parser, ArgParseOption("as", "adapter-seq", "Single adapter sequence as alternative to adapters option.", ARG::STRING));
	addOption(parser, ArgParseOption("ar", "adapter-read-set", "Consider only single read set for adapters.", ARG::STRING));
	addOption(parser, ArgParseOption("ac", "adapter-revcomp", "Consider also reverse complement of each adapter in search."));
	addOption(parser, ArgParseOption("ae", "adapter-trim-end", "Type of removal, see section trim-end modes.", ARG::STRING));
	addOption(parser, ArgParseOption("an", "adapter-tail-length", "Region size for tail trim-end modes. Default: adapter length.", ARG::INTEGER));
	// addOption(parser, ArgParseOption("ah", "adapter-overhang", "Overhang at read ends in right and left modes.", ARG::INTEGER));
	addOption(parser, ArgParseOption("ad", "adapter-relaxed", "Skip restriction to pass read ends in right and left modes."));
	addOption(parser, ArgParseOption("ao", "adapter-min-overlap", "Minimum overlap of adapter and read for removal.", ARG::INTEGER));
	addOption(parser, ArgParseOption("at", "adapter-error-rate", "Error rate threshold for mismatches and gaps.", ARG::DOUBLE));
	addOption(parser, ArgParseOption("am", "adapter-match", "Alignment match score.", ARG::INTEGER));
	addOption(parser, ArgParseOption("ai", "adapter-mismatch", "Alignment mismatch score.", ARG::INTEGER));
	addOption(parser, ArgParseOption("ag", "adapter-gap", "Alignment gap score.", ARG::INTEGER));

	addSection(parser, "Filtering and trimming");
	addOption(parser, ArgParseOption("u", "max-uncalled", "Allowed uncalled bases N for each read.", ARG::INTEGER));
	addOption(parser, ArgParseOption("x", "pre-trim-left", "Trim given number of bases on 5' read end before detection.", ARG::INTEGER));
	addOption(parser, ArgParseOption("y", "pre-trim-right", "Trim specified number of bases on 3' end prior to detection.", ARG::INTEGER));
	addOption(parser, ArgParseOption("k", "post-trim-length", "Trim to specified read length from 3' end after removal.", ARG::INTEGER));
	addOption(parser, ArgParseOption("m", "min-read-length", "Minimum read length to remain after removal.", ARG::INTEGER));

	addSection(parser, "Quality-based trimming");
	addOption(parser, ArgParseOption("q",  "qtrim", "Quality-based trimming mode.", ARG::STRING));
	addOption(parser, ArgParseOption("qf", "qtrim-format", "Quality format.", ARG::STRING));
	addOption(parser, ArgParseOption("qt", "qtrim-threshold", "Minimum quality as threshold for trimming.", ARG::INTEGER));
	// addOption(parser, ArgParseOption("qm", "qtrim-win-mean", "Different threshold for min mean quality of window.", ARG::INTEGER));
	addOption(parser, ArgParseOption("qw", "qtrim-win-size", "Region size for sliding window approach.", ARG::INTEGER));
	addOption(parser, ArgParseOption("qa", "qtrim-post-removal", "Perform quality-based trimming after removal steps."));

	addSection(parser, "Output selection");
	addOption(parser, ArgParseOption("f", "fasta-output", "Prefer non-quality format fasta for output."));
	addOption(parser, ArgParseOption("z", "zip-output", "Direct compression of output files.", ARG::STRING));
	addOption(parser, ArgParseOption("1", "stdout-reads", "Write reads to stdout, tagged and interleaved if needed."));
	addOption(parser, ArgParseOption("j", "length-dist", "Generate length distribution for read output files."));
	addOption(parser, ArgParseOption("s", "single-reads", "Write single reads for too short counterparts in pairs."));
	addOption(parser, ArgParseOption("S", "single-reads-paired", "Write paired single reads with N for short counterparts."));

	addSection(parser, "Logging and tagging");
	addOption(parser, ArgParseOption("l", "align-log", "Print chosen read alignments.", ARG::STRING));
        addOption(parser, ArgParseOption("eve", "everything", "Print all valid alignments between query and read. Additonally mark best alignment with b and log difference between best and second best alignment if it exists"));
        addOption(parser, ArgParseOption("pf", "prefix", "Count number of mismatches and gaps in the prefix of the given adapter/barcode sequences.", ARG::INTEGER));
	addOption(parser, ArgParseOption("o", "stdout-log", "Write statistics to console instead of target log file."));
	addOption(parser, ArgParseOption("g", "removal-tags", "Tag reads that are subject to adapter or barcode removal."));
	addOption(parser, ArgParseOption("e", "number-tags", "Replace read tags by ascending number to save space."));
	addOption(parser, ArgParseOption("d", "random-tags", "Capture read sequence at barcode or adapter N positions."));


	hideOption(parser, "version");

	setAdvanced(parser, "barcodes2");
	setAdvanced(parser, "barcode-tail-length");
	setAdvanced(parser, "barcode-keep");
	setAdvanced(parser, "barcode-unassigned");
	// setAdvanced(parser, "barcode-no-mbv");
	setAdvanced(parser, "barcode-match");
	setAdvanced(parser, "barcode-mismatch");
	setAdvanced(parser, "barcode-gap");

	setAdvanced(parser, "adapters2");
	setAdvanced(parser, "adapter-revcomp");
	setAdvanced(parser, "adapter-tail-length");
	// setAdvanced(parser, "adapter-overhang");
	setAdvanced(parser, "adapter-relaxed");
	setAdvanced(parser, "adapter-read-set");
	setAdvanced(parser, "adapter-match");
	setAdvanced(parser, "adapter-mismatch");
	setAdvanced(parser, "adapter-gap");

	setAdvanced(parser, "post-trim-length");
	setAdvanced(parser, "qtrim-win-size");
	setAdvanced(parser, "qtrim-post-removal");

	setAdvanced(parser, "man-help");
	setAdvanced(parser, "bundle");
	setAdvanced(parser, "stdout-reads");
	setAdvanced(parser, "length-dist");
	setAdvanced(parser, "single-reads-paired");
	setAdvanced(parser, "number-tags");
	setAdvanced(parser, "random-tags");


	setCategory(parser, "Trimming");
	// setRequired(parser, "reads");
	// setMinValue(parser, "threads", "1");

	// setValidValues(parser, "target", "fasta fa fastq fq");
	// setValidValues(parser, "reads", "fasta fa fastq fq");
	// setValidValues(parser, "reads2", "fasta fa fastq fq");
	// setValidValues(parser, "barcode-reads", "fasta fa fastq fq");

	// setValidValues(parser, "barcodes", "fasta fa");
	// setValidValues(parser, "barcodes2", "fasta fa");
	// setValidValues(parser, "adapters", "fasta fa");
	// setValidValues(parser, "adapters2", "fasta fa");

	// setValidValues(parser, "adapter-trim-end", "ANY LEFT RIGHT LTAIL RTAIL");
	// setMinValue(parser, "adapter-tail-length", "1");
	// setMinValue(parser, "adapter-min-overlap", "1");
	// setMinValue(parser, "adapter-error-rate",   "0");
	// setMaxValue(parser, "adapter-error-rate",   "1");

	// setValidValues(parser, "barcode-trim-end", "ANY LEFT RIGHT LTAIL RTAIL");
	// setMinValue(parser, "barcode-tail-length", "1");
	// setMinValue(parser, "barcode-min-overlap", "1");
	// setMinValue(parser, "barcode-error-rate",   "0");
	// setMaxValue(parser, "barcode-error-rate",   "1");

	// setMinValue(parser, "max-uncalled",     "0");
	// setMinValue(parser, "pre-trim-left",    "1");
	// setMinValue(parser, "pre-trim-right",   "1");
	// setMinValue(parser, "post-trim-length", "1");
	// setMinValue(parser, "min-read-length",  "1");
	// setMinValue(parser, "qtrim-threshold",  "0");


	setValidValues(parser, "qtrim", "TAIL WIN BWA");
	setValidValues(parser, "qtrim-format", "sanger solexa i1.3 i1.5 i1.8");
	setValidValues(parser, "align-log", "ALL MOD TAB");
	setValidValues(parser, "zip-output", "GZ BZ2");
	setValidValues(parser, "adapter-read-set", "1 2");

	setDefaultValue(parser, "target",          "flexbarOut");
	setDefaultValue(parser, "threads",         "1");
	setDefaultValue(parser, "bundle",          "256");
	setDefaultValue(parser, "max-uncalled",    "0");
	setDefaultValue(parser, "min-read-length", "18");

	setDefaultValue(parser, "barcode-trim-end",   "LTAIL");
	setDefaultValue(parser, "barcode-error-rate", "0.1");
	setDefaultValue(parser, "barcode-match",      "1");
	setDefaultValue(parser, "barcode-mismatch",   "-1");
	setDefaultValue(parser, "barcode-gap",        "-9");

	setDefaultValue(parser, "adapter-trim-end",    "RIGHT");
	setDefaultValue(parser, "adapter-min-overlap", "3");
	setDefaultValue(parser, "adapter-error-rate",  "0.1");
	// setDefaultValue(parser, "adapter-overhang",    "0");
	setDefaultValue(parser, "adapter-match",       "1");
	setDefaultValue(parser, "adapter-mismatch",    "-1");
	setDefaultValue(parser, "adapter-gap",         "-6");

	setDefaultValue(parser, "qtrim-threshold", "20");
	setDefaultValue(parser, "qtrim-win-size",  "5");


	addTextSection(parser, "TRIM-END MODES");
	addText(parser._toolDoc, "\\fBANY:\\fP   longer side of read remains after removal of overlap", false);
	addText(parser._toolDoc, "\\fBLEFT:\\fP  right side remains after removal, align <= read end",  false);
	addText(parser._toolDoc, "\\fBRIGHT:\\fP left part remains after removal, align >= read start", false);
	addText(parser._toolDoc, "\\fBLTAIL:\\fP consider first n bases of reads in alignment",         false);
	addText(parser._toolDoc, "\\fBRTAIL:\\fP use only last n bases, see tail-length options",       false);

	addTextSection(parser, "EXAMPLES");
	addText(parser._toolDoc, "\\fBflexbar\\fP \\fB-r\\fP reads.fq \\fB-t\\fP target \\fB-b\\fP brc.fa \\fB-be\\fP LTAIL \\fB-a\\fP adp.fa", false);
	addText(parser._toolDoc, "\\fBflexbar\\fP \\fB-r\\fP reads.fq.gz \\fB-q\\fP TAIL \\fB-qf\\fP i1.8 \\fB-a\\fP adp.fa \\fB-ao\\fP 5 \\fB-at\\fP 0.4");
}


void printLocalTime(Options &o){
	time_t t_current;
	time(&t_current);
	*o.out << "Local time:            " << asctime(localtime(&t_current)) << "\n";
}


void parseCmdLine(seqan::ArgumentParser &parser, std::string version, int argc, char const ** argv){

	using namespace std;

	using seqan::ArgumentParser;

	bool useLogFile = true;

	for(int i = 0; i < argc; i++){
		if(strncmp(argv[i], "-o", 2) == 0 || strncmp(argv[i], "--stdout-log",   12) == 0)
			useLogFile = false;
	}
	for(int i = 0; i < argc; i++){
		if(strncmp(argv[i], "-1", 2) == 0 || strncmp(argv[i], "--stdout-reads", 14) == 0)
			useLogFile = true;
	}
	for(int i = 0; i < argc; i++){
		if(strncmp(argv[i], "-h",           2) == 0 ||
		   strncmp(argv[i], "--help",       6) == 0 ||
		   strncmp(argv[i], "--full-help", 11) == 0 ||
		   strncmp(argv[i], "--version",    9) == 0 )
			useLogFile = false;
	}
	if(! useLogFile) cout << endl;


	ArgumentParser::ParseResult res = parse(parser, argc, argv);

	if(res != ArgumentParser::PARSE_OK){

		if(! isSet(parser, "version")){
			cout << endl << getFlexbarURL() << endl;

			if(isSet(parser, "help")){
				cout << "Show advanced options: flexbar -hh\n" << endl;
			}
		}
		else cout << endl;

		exit(res == ArgumentParser::PARSE_ERROR);
	}


	if(isSet(parser, "versions")){
		cout << endl;
		printVersion(parser, cout);
		cout << endl;
		exit(0);
	}
	if(isSet(parser, "cite")){
		cout << endl;
		cout << getFlexbarBanner(version) << endl;
		cout << getFlexbarCitation()      << endl;
		cout << getFlexbarURL()           << endl;
		exit(0);
	}
	if(isSet(parser, "man-help")){
		printHelp(parser, cout, "man", true);
		cout << endl;
		exit(0);
	}
	if(! isSet(parser, "reads")){
		cout << endl;
		printShortHelp(parser);
		cout << endl << getFlexbarURL();
		cerr << "\nPlease specify reads input file.\n" << endl;
		exit(1);
	}
}


void initOptions(Options &o, seqan::ArgumentParser &parser){

	using namespace std;

	bool stdOutReads = isSet(parser, "stdout-reads");
	bool stdOutLog   = isSet(parser, "stdout-log");

	if(stdOutReads) o.useStdout = true;

	if(stdOutLog && ! stdOutReads){
		o.logStdout = true;
		o.out       = &cout;
	}
	else{
		string s;
		getOptionValue(s, parser, "target");
		openOutputFile(o.fstrmOut, s + ".log");

		o.out = &o.fstrmOut;
		*o.out << endl;
	}

	getOptionValue(o.readsFile, parser, "reads");
	checkInputType(o.readsFile, o.format, true);
}


void loadOptions(Options &o, seqan::ArgumentParser &parser){

	using namespace std;
	using namespace flexbar;

	ostream *out = o.out;

	*out << getFlexbarBanner(getVersion(parser)) << endl;
	*out << getFlexbarURL() << endl << endl;

	printLocalTime(o);


	// basic options

	getOptionValue(o.nThreads, parser, "threads");
	*out << "Number of threads:     " << o.nThreads << endl;

	getOptionValue(o.bundleSize, parser, "bundle");
	*out << "Bundled fragments:     " << o.bundleSize << endl << endl;

	getOptionValue(o.targetName, parser, "target");
	*out << "Target name:           " << o.targetName << endl;

	*out << "File type:             ";
	     if(o.format == FASTA)   *out << "fasta";
	else if(o.format == FASTQ)   *out << "fastq";
	*out << endl;


	getOptionValue(o.readsFile, parser, "reads");
	*out << "Reads file:            ";

	if(o.readsFile == "-"){
		*out << "stdin" << endl;
		o.useStdin = true;
	}
	else *out << o.readsFile << endl;

	o.runType = SINGLE;

	if(isSet(parser, "reads2")){
		getOptionValue(o.readsFile2, parser, "reads2");
		*out << "Reads file 2:          " << o.readsFile2 << "   (paired run)" << endl;
		o.runType  = PAIRED;
		o.isPaired = true;

		flexbar::FileFormat fformat;
		checkInputType(o.readsFile2, fformat, false);

		if(o.format != fformat){
			cerr << "\n\n" << "First and second reads file do not have same format.\n" << endl;
			exit(1);
		}
	}


	// barcode and adapter file options

	if(isSet(parser, "barcodes")){

		if(isSet(parser, "barcode-reads")){
			getOptionValue(o.barReadsFile, parser, "barcode-reads");
			*out << "Barcode reads file:    " << o.barReadsFile << endl;

			flexbar::FileFormat fformat;
			checkInputType(o.barReadsFile, fformat, false);

			if(o.format != fformat){
				cerr << "\n\n" << "Barcode reads file does not have same format as reads.\n" << endl;
				exit(1);
			}

			o.barDetect = BARCODE_READ;
		}
		else o.barDetect = WITHIN_READ_REMOVAL;

		getOptionValue(o.barcodeFile, parser, "barcodes");
		*out << "Barcode file:          " << o.barcodeFile << endl;

		     if(o.runType == SINGLE) o.runType = SINGLE_BARCODED;
		else if(o.runType == PAIRED) o.runType = PAIRED_BARCODED;

		if(o.barDetect == WITHIN_READ_REMOVAL && isSet(parser, "barcode-keep")){
			o.barDetect = WITHIN_READ;
		}

		if(isSet(parser, "barcodes2") && o.barDetect != BARCODE_READ && o.isPaired){
			getOptionValue(o.barcode2File, parser, "barcodes2");
			*out << "Barcode file 2:        " << o.barcode2File << endl;

			if(o.barDetect == WITHIN_READ_REMOVAL) o.barDetect = WITHIN_READ_REMOVAL2;
			else if(o.barDetect == WITHIN_READ)    o.barDetect = WITHIN_READ2;
		}
	}

	if(isSet(parser, "adapters")){
		getOptionValue(o.adapterFile, parser, "adapters");
		*out << "Adapter file:          " << o.adapterFile << endl;
		o.adapRm = NORMAL;
		o.useAdapterFile = true;
	}
	else if(isSet(parser, "adapter-seq")){
		getOptionValue(o.adapterSeq, parser, "adapter-seq");
		o.adapRm = NORMAL;
	}

	if(isSet(parser, "adapters2") && o.adapRm == NORMAL && o.isPaired){
		getOptionValue(o.adapter2File, parser, "adapters2");
		*out << "Adapter file 2:        " << o.adapter2File << endl;
		o.adapRm = NORMAL2;
	}
	*out << endl;


	// filtering and trimming options

	getOptionValue(o.maxUncalled, parser, "max-uncalled");
	*out << "max-uncalled:          " << o.maxUncalled << endl;

	if(isSet(parser, "pre-trim-left")){
		getOptionValue(o.cutLen_begin, parser, "pre-trim-left");
		*out << "pre-trim-left:         " << o.cutLen_begin << endl;
	}

	if(isSet(parser, "pre-trim-right")){
		getOptionValue(o.cutLen_end, parser, "pre-trim-right");
		*out << "pre-trim-right:        " << o.cutLen_end << endl;
	}

	if(isSet(parser, "post-trim-length")){
		getOptionValue(o.cutLen_read, parser, "post-trim-length");
		*out << "post-trim-length:      " << o.cutLen_read << endl;
	}

	getOptionValue(o.min_readLen, parser, "min-read-length");
	*out << "min-read-length:       " << o.min_readLen << endl;

	if(o.min_readLen < 1){
		cerr << "\n" << "Minimum read length should be 1 or higher.\n" << endl;
		exit(1);
	}


	// quality-based trimming

	if(isSet(parser, "qtrim") && o.format == FASTQ){

		string qt;
		getOptionValue(qt, parser, "qtrim");

		     if(qt == "TAIL") o.qTrim = TAIL;
 		else if(qt == "WIN")  o.qTrim = WIN;
 		else if(qt == "BWA")  o.qTrim = BWA;

		*out << "qtrim:                 " << qt << endl;

		if(isSet(parser, "qtrim-format")){

			string quality;
			getOptionValue(quality, parser, "qtrim-format");

			     if(quality == "sanger") o.qual = SANGER;
			else if(quality == "solexa") o.qual = SOLEXA;
			else if(quality == "i1.3")   o.qual = ILLUMINA;
			else if(quality == "i1.5")   o.qual = ILLUMINA;
			else if(quality == "i1.8")   o.qual = SANGER;

			*out << "qtrim-format:          " << quality << endl;
		}
		else{
			cerr << "\n\n" << "Specify qtrim-format for quality-based trimming.\n" << endl;
			exit(1);
		}

		getOptionValue(o.qtrimThresh, parser, "qtrim-threshold");

		if(o.qtrimThresh > 0){
			*out << "qtrim-threshold:       " << o.qtrimThresh;

			switch(o.qual){
				case SANGER:      o.qtrimThresh += 33;
					break;
				case SOLEXA:      o.qtrimThresh += 59;
					break;
				case ILLUMINA:  o.qtrimThresh += 64;
			}
			*out << "  (" << o.qtrimThresh << ")" << endl;
		}

		if(o.qTrim == WIN || o.qTrim == WINTAIL){

			// if(isSet(parser, "qtrim-win-mean")){
			// 	getOptionValue(o.qtrimWinMean, parser, "qtrim-win-mean");
			// 	*out << "qtrim-win-mean:        " << o.qtrimWinMean << endl;
			// }

			getOptionValue(o.qtrimWinSize, parser, "qtrim-win-size");
			*out << "qtrim-win-size:        " << o.qtrimWinSize << endl;
		}

		if(isSet(parser, "qtrim-post-removal")) o.qtrimPostRm = true;
	}


	// output, logging and tagging options

	if(isSet(parser, "align-log")){
		getOptionValue(o.logAlignStr, parser, "align-log");

		     if(o.logAlignStr == "ALL") o.logAlign = ALL;
		else if(o.logAlignStr == "TAB") o.logAlign = TAB;
		else if(o.logAlignStr == "MOD") o.logAlign = MOD;
	}

	getOptionValue(o.prefix, parser, "prefix");
	if(isSet(parser, "everything")) o.logEverything = true;

	if(isSet(parser, "zip-output")){
		getOptionValue(o.outCompression, parser, "zip-output");

		if(o.outCompression == "GZ"){
			o.cmprsType = GZ;
			o.outCompression = ".gz";
		}
		else if(o.outCompression == "BZ2"){
			o.cmprsType = BZ2;
			o.outCompression = ".bz2";
		}
	}

	if(isSet(parser, "single-reads")) o.writeSingleReads = true;

	if(isSet(parser, "single-reads-paired")){
		o.writeSingleReadsP = true;
		o.writeSingleReads  = false;
	}

	if(isSet(parser, "fasta-output")) o.switch2Fasta    = true;
	if(isSet(parser, "length-dist"))  o.writeLengthDist = true;
	if(isSet(parser, "number-tags"))  o.useNumberTag    = true;
	if(isSet(parser, "removal-tags")) o.useRemovalTag   = true;
	if(isSet(parser, "random-tags"))  o.randTag         = true;

	*out << endl;


	// barcode options

	if(o.barDetect != BOFF){

		string b_trim_end;
		getOptionValue(b_trim_end, parser, "barcode-trim-end");

		     if(b_trim_end == "LEFT")   o.b_end = LEFT;
		else if(b_trim_end == "RIGHT")  o.b_end = RIGHT;
		else if(b_trim_end == "ANY")    o.b_end = ANY;
		else if(b_trim_end == "LTAIL")  o.b_end = LTAIL;
		else if(b_trim_end == "RTAIL")  o.b_end = RTAIL;
                else if(b_trim_end == "WUN")    o.b_end = WUN;
                else if(b_trim_end == "LTAILS") o.b_end = LTAILS;
		else{
			cerr << "Specified barcode trim-end is unknown!\n" << endl;
			exit(1);
		}
		*out << "barcode-trim-end:      " << b_trim_end << endl;


		if(isSet(parser, "barcode-tail-length")){
			getOptionValue(o.b_tail_len, parser, "barcode-tail-length");
			*out << "barcode-tail-length:   " << o.b_tail_len << endl;
		}

		if(isSet(parser, "barcode-min-overlap")){
			getOptionValue(o.b_min_overlap, parser, "barcode-min-overlap");
			*out << "barcode-min-overlap:   " << o.b_min_overlap << endl;
		}

		getOptionValue(o.b_errorRate, parser, "barcode-error-rate");
		*out << "barcode-error-rate:    " << o.b_errorRate << endl;

		if(o.b_errorRate < 0 || o.b_errorRate >= 1){
			cerr << "\nBarcode error rate should be between 0 and 1.\n" << endl;
			exit(1);
		}

		if(isSet(parser, "barcode-unassigned")) o.writeUnassigned = true;
		// if(isSet(parser, "barcode-no-mbv"))     o.bNoMBV          = true;

		getOptionValue(o.b_match,    parser, "barcode-match");
		getOptionValue(o.b_mismatch, parser, "barcode-mismatch");
		getOptionValue(o.b_gapCost,  parser, "barcode-gap");

		// if(o.bNoMBV){

			*out << "barcode-match:        ";
			if(o.b_match >= 0) *out << " ";
			*out << o.b_match << endl;

			*out << "barcode-mismatch:     ";
			if(o.b_mismatch >= 0) *out << " ";
			*out << o.b_mismatch << endl;

			*out << "barcode-gap:          ";
			if(o.b_gapCost >= 0) *out << " ";
			*out << o.b_gapCost << endl;
		// }
		*out << endl;
	}


	// adapter options

	if(o.adapRm != AOFF){

		string a_trim_end;
		getOptionValue(a_trim_end, parser, "adapter-trim-end");

		if     (a_trim_end == "LEFT")        o.end = LEFT;
		else if(a_trim_end == "RIGHT")       o.end = RIGHT;
		else if(a_trim_end == "ANY")         o.end = ANY;
		else if(a_trim_end == "LTAIL")   o.end = LTAIL;
		else if(a_trim_end == "RTAIL")  o.end = RTAIL;
                else if(a_trim_end == "WUN")    o.end = WUN;
                else if(a_trim_end == "LTAILS") o.end = LTAILS;
		else {
			cerr << "Specified adapter trim-end is unknown!\n" << endl;
			exit(1);
		}
		*out << "adapter-trim-end:      " << a_trim_end << endl;


		if(isSet(parser, "adapter-tail-length")){
			getOptionValue(o.a_tail_len, parser, "adapter-tail-length");
			*out << "adapter-tail-length:   " << o.a_tail_len << endl;
		}

		if(isSet(parser, "adapter-revcomp")){
			*out << "adapter-revcomp:       yes" << endl;
			o.revCompAdapter = true;
		}

		if(isSet(parser, "adapter-relaxed")){
			*out << "adapter-relaxed:       yes" << endl;
			o.relaxRegion = true;
		}

		if(isSet(parser, "adapter-read-set") && o.isPaired && o.adapRm != NORMAL2){
			string a_read_set;
			getOptionValue(a_read_set, parser, "adapter-read-set");
			*out << "adapter-read-set:      " << a_read_set << endl;

			     if(a_read_set == "1") o.adapRm = AONE;
			else if(a_read_set == "2") o.adapRm = ATWO;
		}

		getOptionValue(o.a_min_overlap, parser, "adapter-min-overlap");
		*out << "adapter-min-overlap:   " << o.a_min_overlap << endl;

		getOptionValue(o.a_errorRate, parser, "adapter-error-rate");
		*out << "adapter-error-rate:    " << o.a_errorRate << endl;

		if(o.a_errorRate < 0 || o.a_errorRate >= 1){
			cerr << "\nAdapter error rate should be between 0 and 1.\n" << endl;
			exit(1);
		}

		// getOptionValue(o.a_overhang, parser, "adapter-overhang");
		// *out << "adapter-overhang:      " << o.a_overhang << endl;


		getOptionValue(o.match,    parser, "adapter-match");
		getOptionValue(o.mismatch, parser, "adapter-mismatch");
		getOptionValue(o.gapCost,  parser, "adapter-gap");

		*out << "adapter-match:        ";
		if(o.match >= 0) *out << " ";
		*out << o.match << endl;

		*out << "adapter-mismatch:     ";
		if(o.mismatch >= 0) *out << " ";
		*out << o.mismatch << endl;

		*out << "adapter-gap:          ";
		if(o.gapCost >= 0) *out << " ";
		*out << o.gapCost << "\n" << endl;
	}


	// option compatibility tests

	if(o.cutLen_read != 0 && o.cutLen_read < o.min_readLen){
		o.cutLen_read = 0;
		cerr << "\nOption post-trim-length omitted, as it is shorter than min read length.\n" << endl;
	}

}


#endif
