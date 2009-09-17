//$Id: SprBaggerDecisionTreeApp.cc,v 1.18 2008-04-02 23:36:45 narsky Exp $

#include "StatPatternRecognition/SprExperiment.hh"
#include "StatPatternRecognition/SprDefs.hh"
#include "StatPatternRecognition/SprAbsFilter.hh"
#include "StatPatternRecognition/SprUtils.hh"
#include "StatPatternRecognition/SprBagger.hh"
#include "StatPatternRecognition/SprArcE4.hh"
#include "StatPatternRecognition/SprTrainedBagger.hh"
#include "StatPatternRecognition/SprTopdownTree.hh"
#include "StatPatternRecognition/SprTrainedTopdownTree.hh"
#include "StatPatternRecognition/SprDecisionTree.hh"
#include "StatPatternRecognition/SprTrainedDecisionTree.hh"
#include "StatPatternRecognition/SprEmptyFilter.hh"
#include "StatPatternRecognition/SprAbsReader.hh"
#include "StatPatternRecognition/SprAbsWriter.hh"
#include "StatPatternRecognition/SprDataFeeder.hh"
#include "StatPatternRecognition/SprRWFactory.hh"
#include "StatPatternRecognition/SprTwoClassSignalSignif.hh"
#include "StatPatternRecognition/SprTwoClassIDFraction.hh"
#include "StatPatternRecognition/SprTwoClassTaggerEff.hh"
#include "StatPatternRecognition/SprTwoClassPurity.hh"
#include "StatPatternRecognition/SprTwoClassGiniIndex.hh"
#include "StatPatternRecognition/SprTwoClassCrossEntropy.hh"
#include "StatPatternRecognition/SprTwoClassUniformPriorUL90.hh"
#include "StatPatternRecognition/SprTwoClassBKDiscovery.hh"
#include "StatPatternRecognition/SprTwoClassPunzi.hh"
#include "StatPatternRecognition/SprIntegerBootstrap.hh"
#include "StatPatternRecognition/SprStringParser.hh"
#include "StatPatternRecognition/SprCrossValidator.hh"
#include "StatPatternRecognition/SprAverageLoss.hh"
#include "StatPatternRecognition/SprLoss.hh"
#include "StatPatternRecognition/SprClass.hh"
#include "StatPatternRecognition/SprClassifierReader.hh"
#include "StatPatternRecognition/SprAbsVarTransformer.hh"
#include "StatPatternRecognition/SprVarTransformerReader.hh"
#include "StatPatternRecognition/SprTransformerFilter.hh"

#include <stdlib.h>
#include <unistd.h>
#include <iostream>
#include <vector>
#include <set>
#include <string>
#include <memory>
#include <iomanip>

using namespace std;


void help(const char* prog) 
{
  cout << "Usage:  " << prog 
       << " training_data_file" << endl;
  cout << "\t Options: " << endl;
  cout << "\t-h --- help                                        " << endl;
  cout << "\t-j use regular tree instead of faster topdown tree " << endl;
  cout << "\t-k discrete decision tree output (default=continuous)"<< endl;
  cout << "\t-o output Tuple file                               " << endl;
  cout << "\t-a input ascii file mode (see SprSimpleReader.hh)  " << endl;
  cout << "\t-A save output data in ascii instead of Root       " << endl;
  cout << "\t-n number of Bagger training cycles                " << endl;
  cout << "\t-l minimal number of entries per tree leaf (def=1) " << endl;
  cout << "\t-s max number of sampled features (def=0 no sampling)"<< endl;
  cout << "\t-y list of input classes (see SprAbsFilter.hh)     " << endl;
  cout << "\t-Q apply variable transformation saved in file     " << endl;
  cout << "\t-B use a version of Breiman's arc-x4 algorithm     " << endl;
  cout << "\t-v verbose level (0=silent default,1,2)            " << endl;
  cout << "\t-f store trained Bagger to file                    " << endl;
  cout << "\t-F generate code for AdaBoost and store to file    " << endl;
  cout << "\t-c criterion for decision tree optimization (def=5)" << endl;
  cout << "\t-C criterion to be monitored on test data (def=1)  " << endl;
  cout << "\t\t 1 = correctly classified fraction               " << endl;
  cout << "\t\t 2 = signal significance s/sqrt(s+b)             " << endl;
  cout << "\t\t 3 = purity s/(s+b)                              " << endl;
  cout << "\t\t 4 = tagger efficiency Q                         " << endl;
  cout << "\t\t 5 = Gini index                                  " << endl;
  cout << "\t\t 6 = cross-entropy                               " << endl;
  cout << "\t\t 7 = 90% Bayesian upper limit with uniform prior " << endl;
  cout << "\t\t 8 = discovery potential 2*(sqrt(s+b)-sqrt(b))   " << endl;
  cout << "\t\t 9 = Punzi's sensitivity s/(0.5*nSigma+sqrt(b))  " << endl;
  cout << "\t\t -P background normalization factor for Punzi FOM" << endl;
  cout << "\t-g per-event loss for (cross-)validation           " << endl;
  cout << "\t\t 1 - quadratic loss (y-f(x))^2                   " << endl;
  cout << "\t\t 2 - exponential loss exp(-y*f(x))               " << endl;
  cout << "\t\t 3 - misid fraction                              " << endl;
  cout << "\t-i count splits on input variables                 " << endl;
  cout << "\t-r resume training for Bagger stored in file       " << endl;
  cout << "\t-N use only N trees from trained Bagger            " << endl;
  cout << "\t-K keep this fraction in training set and          " << endl;
  cout << "\t\t put the rest into validation set                " << endl;
  cout << "\t-D randomize training set split-up                 " << endl;
  cout << "\t-G generate seed from time of day for bootstrap    " << endl;
  cout << "\t-t read validation/test data from a file           " << endl;
  cout << "\t\t (must be in same format as input data!!!        " << endl;
  cout << "\t-d frequency of print-outs for validation data     " << endl;
  cout << "\t-w scale all signal weights by this factor         " << endl;
  cout << "\t-V include only these input variables              " << endl;
  cout << "\t-z exclude input variables from the list           " << endl;
  cout << "\t-Z exclude input variables from the list, "
       << "but put them in the output file " << endl;
  cout << "\t\t Variables must be listed in quotes and separated by commas." 
       << endl;
  cout << "\t-x cross-validate by splitting data into a given "
       << "number of pieces" << endl;
  cout << "\t-q a set of minimal node sizes for cross-validation" << endl;
  cout << "\t\t Node sizes must be listed in quotes and separated by commas." 
       << endl;
}


int main(int argc, char ** argv)
{
  // check command line
  if( argc < 2 ) {
    help(argv[0]);
    return 1;
  }

  // init
  string tupleFile;
  int readMode = 0;
  SprRWFactory::DataType writeMode = SprRWFactory::Root;
  unsigned cycles = 0;
  unsigned nmin = 1;
  int verbose = 0;
  string outFile;
  string codeFile;
  string resumeFile;
  int iCritOpt = 5;
  int iCritMon = 1;
  string valFile;
  unsigned valPrint = 0;
  bool scaleWeights = false;
  double sW = 1.;
  int nFeaturesToSample = 0;
  bool countTreeSplits = false;
  string includeList, excludeList;
  unsigned nCross = 0;
  string nodeValidationString;
  bool useTopdown = true;
  bool discrete = false;
  int iLoss = 0;
  string inputClassesString;
  int baggerMode = 1; // 1 = standard  2 = ArcX4
  double bW = 1.;
  string stringVarsDoNotFeed;
  bool split = false;
  double splitFactor = 0;
  bool splitRandomize = false;
  bool initBootstrapFromTimeOfDay = false;
  string transformerFile;
  unsigned nUsedTrees = 0;

  // decode command line
  int c;
  extern char* optarg;
  //  extern int optind;
  while( (c = getopt(argc,argv,"hjko:a:An:l:s:y:Q:Bv:f:F:c:C:P:g:ir:N:K:DGt:d:w:V:z:Z:x:q:")) 
	 != EOF ) {
    switch( c )
      {
      case 'h' :
	help(argv[0]);
	return 1;
      case 'j' :
	useTopdown = false;
	break;
      case 'k' :
	discrete = true;
	break;
      case 'o' :
	tupleFile = optarg;
	break;
      case 'a' :
	readMode = (optarg==0 ? 0 : atoi(optarg));
	break;
      case 'A' :
	writeMode = SprRWFactory::Ascii;
	break;
      case 'n' :
	cycles = (optarg==0 ? 1 : atoi(optarg));
	break;
      case 'l' :
	nmin = (optarg==0 ? 1 : atoi(optarg));
	break;
      case 's' :
	nFeaturesToSample = (optarg==0 ? 0 : atoi(optarg));
	break;
      case 'y' :
	inputClassesString = optarg;
	break;
      case 'Q' :
        transformerFile = optarg;
        break;
      case 'B' :
	baggerMode = 2;
	break;
      case 'v' :
	verbose = (optarg==0 ? 0 : atoi(optarg));
	break;
      case 'f' :
	outFile = optarg;
	break;
      case 'F' :
	codeFile = optarg;
	break;
      case 'c' :
        iCritOpt = (optarg==0 ? 5 : atoi(optarg));
        break;
      case 'C' :
        iCritMon = (optarg==0 ? 1 : atoi(optarg));
        break;
      case 'P' :
	bW = (optarg==0 ? 1 : atof(optarg));
	break;
      case 'g' :
        iLoss = (optarg==0 ? 0 : atoi(optarg));
        break;
      case 'i' :
	countTreeSplits = true;
	break;
      case 'r' :
	resumeFile = optarg;
	break;
      case 'N' :
	nUsedTrees = (optarg==0 ? 0 : atoi(optarg));
	break;
      case 'K' :
        split = true;
        splitFactor = (optarg==0 ? 0 : atof(optarg));
        break;
      case 'D' :
        splitRandomize = true;
        break;
      case 'G' :
	initBootstrapFromTimeOfDay = true;
	break;
      case 't' :
	valFile = optarg;
	break;
      case 'd' :
	valPrint = (optarg==0 ? 0 : atoi(optarg));
	break;
      case 'w' :
	if( optarg != 0 ) {
	  scaleWeights = true;
	  sW = atof(optarg);
	}
	break;
      case 'V' :
	includeList = optarg;
	break;
      case 'z' :
	excludeList = optarg;
	break;
      case 'Z' :
	stringVarsDoNotFeed = optarg;
	break;
      case 'x' :
	nCross = (optarg==0 ? 0 : atoi(optarg));
	break;
      case 'q' :
	nodeValidationString = optarg;
	break;
      }
  }

  // There has to be 1 argument after all options.
  string trFile = argv[argc-1];
  if( trFile.empty() ) {
    cerr << "No training file is specified." << endl;
    return 1;
  }

  // make reader
  SprRWFactory::DataType inputType 
    = ( readMode==0 ? SprRWFactory::Root : SprRWFactory::Ascii );
  auto_ptr<SprAbsReader> reader(SprRWFactory::makeReader(inputType,readMode));

  // include variables
  set<string> includeSet;
  if( !includeList.empty() ) {
    vector<vector<string> > includeVars;
    SprStringParser::parseToStrings(includeList.c_str(),includeVars);
    assert( !includeVars.empty() );
    for( int i=0;i<includeVars[0].size();i++ ) 
      includeSet.insert(includeVars[0][i]);
    if( !reader->chooseVars(includeSet) ) {
      cerr << "Unable to include variables in training set." << endl;
      return 2;
    }
    else {
      cout << "Following variables have been included in optimization: ";
      for( set<string>::const_iterator 
	     i=includeSet.begin();i!=includeSet.end();i++ )
	cout << "\"" << *i << "\"" << " ";
      cout << endl;
    }
  }

  // exclude variables
  set<string> excludeSet;
  if( !excludeList.empty() ) {
    vector<vector<string> > excludeVars;
    SprStringParser::parseToStrings(excludeList.c_str(),excludeVars);
    assert( !excludeVars.empty() );
    for( int i=0;i<excludeVars[0].size();i++ ) 
      excludeSet.insert(excludeVars[0][i]);
    if( !reader->chooseAllBut(excludeSet) ) {
      cerr << "Unable to exclude variables from training set." << endl;
      return 2;
    }
    else {
      cout << "Following variables have been excluded from optimization: ";
      for( set<string>::const_iterator 
	     i=excludeSet.begin();i!=excludeSet.end();i++ )
	cout << "\"" << *i << "\"" << " ";
      cout << endl;
    }
  }

  // read training data from file
  auto_ptr<SprAbsFilter> filter(reader->read(trFile.c_str()));
  if( filter.get() == 0 ) {
    cerr << "Unable to read data from file " << trFile.c_str() << endl;
    return 2;
  }
  vector<string> vars;
  filter->vars(vars);
  cout << "Read data from file " << trFile.c_str() 
       << " for variables";
  for( int i=0;i<vars.size();i++ ) 
    cout << " \"" << vars[i].c_str() << "\"";
  cout << endl;
  cout << "Total number of points read: " << filter->size() << endl;

  // filter training data by class
  vector<SprClass> inputClasses;
  if( !filter->filterByClass(inputClassesString.c_str()) ) {
    cerr << "Cannot choose input classes for string " 
	 << inputClassesString << endl;
    return 2;
  }
  filter->classes(inputClasses);
  assert( inputClasses.size() > 1 );
  cout << "Training data filtered by class." << endl;
  for( int i=0;i<inputClasses.size();i++ ) {
    cout << "Points in class " << inputClasses[i] << ":   " 
	 << filter->ptsInClass(inputClasses[i]) << endl;
  }

  // scale weights
  if( scaleWeights ) {
    cout << "Signal weights are multiplied by " << sW << endl;
    filter->scaleWeights(inputClasses[1],sW);
  }

  // read validation data from file
  auto_ptr<SprAbsFilter> valFilter;
  if( split && !valFile.empty() ) {
    cerr << "Unable to split training data and use validation data "
         << "from a separate file." << endl;
    return 2;
  }
  if( split && valPrint!=0 ) {
    cout << "Splitting training data with factor " << splitFactor << endl;
    if( splitRandomize )
      cout << "Will use randomized splitting." << endl;
    vector<double> weights;
    SprData* splitted = filter->split(splitFactor,weights,splitRandomize);
    if( splitted == 0 ) {
      cerr << "Unable to split training data." << endl;
      return 2;
    }
    bool ownData = true;
    valFilter.reset(new SprEmptyFilter(splitted,weights,ownData));
    cout << "Training data re-filtered:" << endl;
    for( int i=0;i<inputClasses.size();i++ ) {
      cout << "Points in class " << inputClasses[i] << ":   "
           << filter->ptsInClass(inputClasses[i]) << endl;
    }
  }
  if( !valFile.empty() && valPrint!=0 ) {
    auto_ptr<SprAbsReader> 
      valReader(SprRWFactory::makeReader(inputType,readMode));
    if( !includeSet.empty() ) {
      if( !valReader->chooseVars(includeSet) ) {
	cerr << "Unable to include variables in validation set." << endl;
	return 2;
      }
    }
    if( !excludeSet.empty() ) {
      if( !valReader->chooseAllBut(excludeSet) ) {
	cerr << "Unable to exclude variables from validation set." << endl;
	return 2;
      }
    }
    valFilter.reset(valReader->read(valFile.c_str()));
    if( valFilter.get() == 0 ) {
      cerr << "Unable to read data from file " << valFile.c_str() << endl;
      return 2;
    }
    vector<string> valVars;
    valFilter->vars(valVars);
    cout << "Read validation data from file " << valFile.c_str() 
	 << " for variables";
    for( int i=0;i<valVars.size();i++ ) 
      cout << " \"" << valVars[i].c_str() << "\"";
    cout << endl;
    cout << "Total number of points read: " << valFilter->size() << endl;
  }

  // filter validation data by class
  if( valFilter.get() != 0 ) {
    if( !valFilter->filterByClass(inputClassesString.c_str()) ) {
      cerr << "Cannot choose input classes for string " 
	   << inputClassesString << endl;
      return 2;
    }
    valFilter->classes(inputClasses);
    cout << "Validation data filtered by class." << endl;
    for( int i=0;i<inputClasses.size();i++ ) {
      cout << "Points in class " << inputClasses[i] << ":   " 
	   << valFilter->ptsInClass(inputClasses[i]) << endl;
    }
  }

  // scale weights
  if( scaleWeights && valFilter.get()!=0 )
    valFilter->scaleWeights(inputClasses[1],sW);

  // apply transformation of variables to training and test data
  auto_ptr<SprAbsFilter> garbage_train, garbage_valid;
  if( !transformerFile.empty() ) {
    const SprAbsVarTransformer* t 
      = SprVarTransformerReader::read(transformerFile.c_str());
    if( t == 0 ) {
      cerr << "Unable to read VarTransformer from file "
           << transformerFile.c_str() << endl;
      return 2;
    }
    SprTransformerFilter* t_train = new SprTransformerFilter(filter.get());
    SprTransformerFilter* t_valid = 0;
    if( valFilter.get() != 0 )
      t_valid = new SprTransformerFilter(valFilter.get());
    bool replaceOriginalData = true;
    if( !t_train->transform(t,replaceOriginalData) ) {
      cerr << "Unable to apply VarTransformer to training data." << endl;
      return 2;
    }
    if( t_valid!=0 && !t_valid->transform(t,replaceOriginalData) ) {
      cerr << "Unable to apply VarTransformer to validation data." << endl;
      return 2;
    }
    cout << "Variable transformation from file "
         << transformerFile.c_str() << " has been applied to "
         << "training and validation data." << endl;
    garbage_train.reset(filter.release());
    garbage_valid.reset(valFilter.release());
    filter.reset(t_train);
    valFilter.reset(t_valid);
  }

  // make optimization criterion
  auto_ptr<SprAbsTwoClassCriterion> critOpt;
  switch( iCritOpt )
    {
    case 1 :
      critOpt.reset(new SprTwoClassIDFraction);
      cout << "Optimization criterion set to "
           << "Fraction of correctly classified events " << endl;
      break;
    case 2 :
      critOpt.reset(new SprTwoClassSignalSignif);
      cout << "Optimization criterion set to "
           << "Signal significance S/sqrt(S+B) " << endl;
      break;
    case 3 :
      critOpt.reset(new SprTwoClassPurity);
      cout << "Optimization criterion set to "
           << "Purity S/(S+B) " << endl;
      break;
    case 4 :
      critOpt.reset(new SprTwoClassTaggerEff);
      cout << "Optimization criterion set to "
           << "Tagging efficiency Q = e*(1-2w)^2 " << endl;
      break;
    case 5 :
      critOpt.reset(new SprTwoClassGiniIndex);
      cout << "Optimization criterion set to "
	   << "Gini index  -1+p^2+q^2 " << endl;
      break;
    case 6 :
      critOpt.reset(new SprTwoClassCrossEntropy);
      cout << "Optimization criterion set to "
	   << "Cross-entropy p*log(p)+q*log(q) " << endl;
      break;
    case 7 :
      critOpt.reset(new SprTwoClassUniformPriorUL90);
      cout << "Optimization criterion set to "
           << "Inverse of 90% Bayesian upper limit with uniform prior" << endl;
      break;
    case 8 :
      critOpt.reset(new SprTwoClassBKDiscovery);
      cout << "Optimization criterion set to "
	   << "Discovery potential 2*(sqrt(S+B)-sqrt(B))" << endl;
      break;
    case 9 :
      critOpt.reset(new SprTwoClassPunzi(bW));
      cout << "Optimization criterion set to "
	   << "Punzi's sensitivity S/(0.5*nSigma+sqrt(B))" << endl;
      break;
    default :
      cerr << "Unable to make initialization criterion." << endl;
      return 3;
    }

  // check criterion vs classifier
  if( baggerMode==2 && !critOpt->symmetric() ) {
    cerr << "Unable to use arc-e4 with an asymmetric criterion." << endl;
    return 3;
  }

  // make monitoring criterion
  auto_ptr<SprAbsTwoClassCriterion> critMon;
  switch( iCritMon )
    {
    case 1 :
      critMon.reset(new SprTwoClassIDFraction);
      cout << "Monitoring criterion set to "
           << "Fraction of correctly classified events " << endl;
      break;
    case 2 :
      critMon.reset(new SprTwoClassSignalSignif);
      cout << "Monitoring criterion set to "
           << "Signal significance S/sqrt(S+B) " << endl;
      break;
    case 3 :
      critMon.reset(new SprTwoClassPurity);
      cout << "Monitoring criterion set to "
           << "Purity S/(S+B) " << endl;
      break;
    case 4 :
      critMon.reset(new SprTwoClassTaggerEff);
      cout << "Monitoring criterion set to "
           << "Tagging efficiency Q = e*(1-2w)^2 " << endl;
      break;
    case 5 :
      critMon.reset(new SprTwoClassGiniIndex);
      cout << "Monitoring criterion set to "
	   << "Gini index  -1+p^2+q^2 " << endl;
      break;
    case 6 :
      critMon.reset(new SprTwoClassCrossEntropy);
      cout << "Monitoring criterion set to "
	   << "Cross-entropy p*log(p)+q*log(q) " << endl;
      break;
    case 7 :
      critMon.reset(new SprTwoClassUniformPriorUL90);
      cout << "Monitoring criterion set to "
           << "Inverse of 90% Bayesian upper limit with uniform prior" << endl;
      break;
    case 8 :
      critMon.reset(new SprTwoClassBKDiscovery);
      cout << "Monitoring criterion set to "
	   << "Discovery potential 2*(sqrt(S+B)-sqrt(B))" << endl;
      break;
    case 9 :
      critMon.reset(new SprTwoClassPunzi(bW));
      cout << "Monitoring criterion set to "
	   << "Punzi's sensitivity S/(0.5*nSigma+sqrt(B))" << endl;
      break;
    default :
      cerr << "Unable to make initialization criterion." << endl;
      return 3;
    }

  // make per-event loss
  auto_ptr<SprAverageLoss> loss;
  switch( iLoss )
    {
    case 1 :
      loss.reset(new SprAverageLoss(&SprLoss::quadratic));
      cout << "Per-event loss set to "
           << "Quadratic loss (y-f(x))^2 " << endl;
      break;
    case 2 :
      loss.reset(new SprAverageLoss(&SprLoss::purity_ratio));
      cout << "Per-event loss set to "
           << "Exponential loss exp(-y*f(x)) " << endl;
      break;
    case 3 :
      loss.reset(new SprAverageLoss(&SprLoss::correct_id));
      cout << "Per-event loss set to "
	   << "Misid rate int(y==f(x)) " << endl;
      break;
    default :
      cout << "No per-event loss is chosen. Will use the default." << endl;
      break;
    }

  // make bootstrap for resampling input features
  auto_ptr<SprIntegerBootstrap> bootstrap;
  if( nFeaturesToSample > filter->dim() ) 
    nFeaturesToSample = filter->dim();
  if( nFeaturesToSample > 0 ) {
    bootstrap.reset(new SprIntegerBootstrap(filter->dim(),nFeaturesToSample));
    if( !resumeFile.empty() || initBootstrapFromTimeOfDay ) 
      bootstrap->init(-1);
  }

  // make decision tree
  bool doMerge = !critOpt->symmetric();
  if( doMerge ) useTopdown = false;
  auto_ptr<SprDecisionTree> tree;
  if( useTopdown ) {
    tree.reset(new SprTopdownTree(filter.get(),critOpt.get(),nmin,
				  discrete,bootstrap.get()));
  }
  else {
    tree.reset(new SprDecisionTree(filter.get(),critOpt.get(),nmin,doMerge,
				   discrete,bootstrap.get()));
  }
  if( countTreeSplits ) tree->startSplitCounter();
  tree->useFastSort();

  // if cross-validation requested, cross-validate and exit
  if( nCross > 0 ) {
    // message
    cout << "Will cross-validate by dividing training data into " 
	 << nCross << " subsamples." << endl;
    vector<vector<int> > nodeMinSize;

    // decode validation string
    if( !nodeValidationString.empty() )
      SprStringParser::parseToInts(nodeValidationString.c_str(),nodeMinSize);
    else {
      nodeMinSize.resize(1);
      nodeMinSize[0].push_back(nmin);
    }
    if( nodeMinSize.empty() || nodeMinSize[0].empty() ) {
      cerr << "Unable to determine node size for cross-validation." << endl;
      return 4;
    }
    else {
      cout << "Will cross-validate for trees with minimal node sizes: ";
      for( int i=0;i<nodeMinSize[0].size();i++ )
	cout << nodeMinSize[0][i] << " ";
      cout << endl;
    }

    // loop over nodes to prepare classifiers
    vector<SprDecisionTree*> trees(nodeMinSize[0].size());
    vector<SprAbsClassifier*> classifiers(nodeMinSize[0].size());
    for( int i=0;i<nodeMinSize[0].size();i++ ) {
      SprDecisionTree* tree1 = 0;
      if( useTopdown ) {
	tree1 = new SprTopdownTree(filter.get(),critOpt.get(),
				   nodeMinSize[0][i],
				   discrete,bootstrap.get());
      }
      else {
	tree1 = new SprDecisionTree(filter.get(),critOpt.get(),
				    nodeMinSize[0][i],
				    doMerge,discrete,bootstrap.get());
      }
      tree1->useFastSort();
      SprBagger* bagger1 = 0;
      if(      baggerMode == 1 )
	bagger1 = new SprBagger(filter.get(),cycles,discrete);
      else if( baggerMode == 2 )
	bagger1 = new SprArcE4(filter.get(),cycles,discrete);
      if( initBootstrapFromTimeOfDay 
	  && !bagger1->initBootstrapFromTimeOfDay() ) {
	cerr << "Unable to generate seed from time of day for Bagger." << endl;
	return 4;
      }
      if( !bagger1->addTrainable(tree1) ) {
	cerr << "Unable to add decision tree to Bagger for CV." << endl;
	for( int j=0;j<trees.size();j++ ) {
	  delete trees[j];
	  delete classifiers[j];
	}
	return 4;
      }
      trees[i] = tree1;
      classifiers[i] = bagger1;
    }

    // cross-validate
    vector<SprValueWithError> cvFom;
    SprCrossValidator cv(filter.get(),nCross);
    if( !cv.validate(critMon.get(),loss.get(),classifiers,
		     inputClasses[0],inputClasses[1],cvFom,false,verbose) ) {
      cerr << "Unable to cross-validate." << endl;
      for( int j=0;j<trees.size();j++ ) {
	delete trees[j];
	delete classifiers[j];
      }
      return 4;
    }
    else {
      cout << "Cross-validated FOMs:" << endl;
      for( int i=0;i<cvFom.size();i++ ) {
	cout << "Node size=" << setw(8) << nodeMinSize[0][i] 
	     << "      FOM=" << setw(10) << cvFom[i].first 
	     << " +- " << setw(10) << cvFom[i].second << endl;
      }
    }

    // cleanup
    for( int j=0;j<trees.size();j++ ) {
      delete trees[j];
      delete classifiers[j];
    }

    // normal exit
    return 0;
  }// end cross-validation

  // make Bagger
  auto_ptr<SprBagger> bagger;
  if(      baggerMode == 1 )
    bagger.reset(new SprBagger(filter.get(),cycles,discrete));
  else if( baggerMode == 2 )
    bagger.reset(new SprArcE4(filter.get(),cycles,discrete));

  // set seed for bootstrap if necessary
  if( initBootstrapFromTimeOfDay && !bagger->initBootstrapFromTimeOfDay() ) {
    cerr << "Unable to generate seed from time of day for Bagger." << endl;
    return 4;
  }

  // set validation
  if( valFilter.get()!=0 && !valFilter->empty() )
    bagger->setValidation(valFilter.get(),valPrint,critMon.get(),loss.get());

  // read saved Bagger from file
  if( !resumeFile.empty() ) {
    if( !SprClassifierReader::readTrainable(resumeFile.c_str(),
                                            bagger.get(),verbose) ) {
      cerr << "Failed to read saved Bagger from file "
           << resumeFile.c_str() << endl;
      return 5;
    }
    cout << "Read saved Bagger from file " << resumeFile.c_str()
	 << " with " << bagger->nTrained() << " trained classifiers." 
	 << endl;
  }

  // add trainable tree
  if( !bagger->addTrainable(tree.get()) ) {
    cerr << "Unable to add decision tree to Bagger." << endl;
    return 6;
  }

  // train
  if( !bagger->train(verbose) )
    cerr << "Bagger terminated with error." << endl;
  if( bagger->nTrained() == 0 ) {
    cerr << "Unable to train Bagger." << endl;
    return 7;
  }
  else {
    cout << "Bagger finished training with " << bagger->nTrained() 
	 << " classifiers." << endl;
  }

  // save trained Bagger
  if( !outFile.empty() ) {
    if( !bagger->store(outFile.c_str()) ) {
      cerr << "Cannot store Bagger in file " << outFile.c_str() << endl;
      return 8;
    }
  }

  // print out counted splits
  if( countTreeSplits ) tree->printSplitCounter(cout);

  // make a trained Bagger
  auto_ptr<SprTrainedBagger> trainedBagger(bagger->makeTrained());
  if( trainedBagger.get() == 0 ) {
    cerr << "Unable to get trained Bagger." << endl;
    return 7;
  }

  // reduce the number of trees if requested by the user
  if( nUsedTrees > 0 ) {
    cout << "Will use " << nUsedTrees << " trees instead of " 
	 << trainedBagger->nClassifiers() << " trained by Bagger." << endl;
    trainedBagger->useNClassifiers(nUsedTrees);
  }

  // store code into file
  if( !codeFile.empty() ) {
    if( !trainedBagger->storeCode(codeFile.c_str()) ) {
      cerr << "Unable to store code for trained Bagger." << endl;
      return 8;
    }
  }

  // make histogram if requested
  if( tupleFile.empty() ) 
    return 0;

  // make a writer
  auto_ptr<SprAbsWriter> tuple(SprRWFactory::makeWriter(writeMode,"training"));
  if( !tuple->init(tupleFile.c_str()) ) {
    cerr << "Unable to open output file " << tupleFile.c_str() << endl;
    return 9;
  }

  // determine if certain variables are to be excluded from usage,
  // but included in the output storage file (-Z option)
  string printVarsDoNotFeed;
  vector<vector<string> > varsDoNotFeed;
  SprStringParser::parseToStrings(stringVarsDoNotFeed.c_str(),varsDoNotFeed);
  vector<unsigned> mapper;
  for( int d=0;d<vars.size();d++ ) {
    if( varsDoNotFeed.empty() ||
        (find(varsDoNotFeed[0].begin(),varsDoNotFeed[0].end(),vars[d])
	 ==varsDoNotFeed[0].end()) ) {
      mapper.push_back(d);
    }
    else {
      printVarsDoNotFeed += ( printVarsDoNotFeed.empty() ? "" : ", " );
      printVarsDoNotFeed += vars[d];
    }
  }
  if( !printVarsDoNotFeed.empty() ) {
    cout << "The following variables are not used in the algorithm, " 
         << "but will be included in the output file: " 
         << printVarsDoNotFeed.c_str() << endl;
  }

  // feed
  SprDataFeeder feeder(filter.get(),tuple.get(),mapper);
  feeder.addClassifier(trainedBagger.get(),"bag");
  if( !feeder.feed(1000) ) {
    cerr << "Cannot feed data into file " << tupleFile.c_str() << endl;
    return 10;
  }

  // exit
  return 0;
}
