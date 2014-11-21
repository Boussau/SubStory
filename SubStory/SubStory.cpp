//
// File: SubStory.cpp
// Created by: Bastien Boussau
// Created on: Nov Fri 21 2014
//

/*
Copyright or © or Copr. Bio++ Development Team

This software is a computer program whose purpose is to estimate
phylogenies and evolutionary parameters from a dataset according to
the maximum likelihood principle.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

// From the STL:
#include <iostream>
#include <iomanip>
#include <map>

using namespace std;

#include <Bpp/App/ApplicationTools.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/Prob/GammaDiscreteDistribution.h>
#include <Bpp/Numeric/DataTable.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Numeric/Stat/StatTools.h>
#include <Bpp/Text/KeyvalTools.h>
#include <Bpp/Text/TextTools.h>

// From SeqLib:
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/AlphabetIndex/UserAlphabetIndex1.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/SequenceContainerTools.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>

// From PhylLib:
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/Likelihood.all>
#include <Bpp/Phyl/Mapping.all>
#include <Bpp/Phyl/PatternTools.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Io/Nhx.h>
#include <Bpp/Phyl/Io/BppOSubstitutionModelFormat.h>
#include <Bpp/Phyl/Model.all>
#include <Bpp/Phyl/Model/MarkovModulatedSubstitutionModel.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Model/SubstitutionModelSet.h>
#include <Bpp/Phyl/Model/SubstitutionModelSetTools.h>



using namespace bpp;


/******************************************************************************/

void help()
{
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
  (*ApplicationTools::message << "substory parameter1_name=parameter1_value ").endLine();
  (*ApplicationTools::message << "      parameter2_name=parameter2_value ... param=option_file").endLine();
  (*ApplicationTools::message).endLine();
  (*ApplicationTools::message << "  Refer to the Bio++ Manual for a list of available options.").endLine();
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
}

int main(int args, char ** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*     Bio++ Substitution History Reconstruction, version 2.2.0    *" << endl;
  cout << "* Authors: B. Boussau                       Created on: 21/11/14 *" << endl;
  cout << "*                                           Last Modif: 21/11/14 *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  if (args == 1)
  {
    help();
    return 0;
  }
  
  try {

  BppApplication substory(args, argv, "substory");
  substory.startTimer();

  Alphabet* alphabet = SequenceApplicationTools::getAlphabet(substory.getParams(), "", false);
  auto_ptr<GeneticCode> gCode;
  CodonAlphabet* codonAlphabet = dynamic_cast<CodonAlphabet*>(alphabet);
  if (codonAlphabet) {
    string codeDesc = ApplicationTools::getStringParameter("genetic_code", substory.getParams(), "Standard", "", true, true);
    ApplicationTools::displayResult("Genetic Code", codeDesc);
      
    gCode.reset(SequenceApplicationTools::getGeneticCode(codonAlphabet->getNucleicAlphabet(), codeDesc));
  }

  VectorSiteContainer* allSites = SequenceApplicationTools::getSiteContainer(alphabet, substory.getParams());
  
  VectorSiteContainer* sites = SequenceApplicationTools::getSitesToAnalyse(* allSites, substory.getParams(), "", true, false);
  delete allSites;

  ApplicationTools::displayResult("Number of sequences", TextTools::toString(sites->getNumberOfSequences()));
  ApplicationTools::displayResult("Number of sites", TextTools::toString(sites->getNumberOfSites()));
  
  // Get the initial tree
  Tree* tree = PhylogeneticsApplicationTools::getTree(substory.getParams());
  ApplicationTools::displayResult("Number of leaves", TextTools::toString(tree->getNumberOfLeaves()));
  
  string treeWIdPath = ApplicationTools::getAFilePath("output.tree_ids.file", substory.getParams(), false, false);
  if (treeWIdPath != "none")
  {
    TreeTemplate<Node> ttree(*tree);
    vector<Node *> nodes = ttree.getNodes();
    for(unsigned int i = 0; i < nodes.size(); i++)
    {
      if(nodes[i]->isLeaf())
        nodes[i]->setName(TextTools::toString(nodes[i]->getId()) + "_" + nodes[i]->getName());
      else
        nodes[i]->setBranchProperty("NodeId", BppString(TextTools::toString(nodes[i]->getId())));
    }
    Newick treeWriter;
    treeWriter.enableExtendedBootstrapProperty("NodeId");
    ApplicationTools::displayResult("Writing tagged tree to", treeWIdPath);
    treeWriter.write(ttree, treeWIdPath);
    delete tree;
    cout << "BppAncestor's done." << endl;
    exit(0);
  }

  bool checkTree = ApplicationTools::getBooleanParameter("input.tree.check_root", substory.getParams(), true, "", true, false);

  DRTreeLikelihood *tl;
  string nhOpt = ApplicationTools::getStringParameter("nonhomogeneous", substory.getParams(), "no", "", true, false);
  ApplicationTools::displayResult("Heterogeneous model", nhOpt);

  SubstitutionModel    *model    = 0;
  SubstitutionModelSet *modelSet = 0;
  DiscreteDistribution *rDist    = 0;
  size_t nbStates;

  if (nhOpt == "no")
  {  
    model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, gCode.get(), sites, substory.getParams());
    if (model->getName() != "RE08") SiteContainerTools::changeGapsToUnknownCharacters(*sites);
    if (model->getNumberOfStates() > model->getAlphabet()->getSize())
    {
      //Markov-modulated Markov model!
      rDist = new ConstantRateDistribution();
    }
    else
    {
      rDist = PhylogeneticsApplicationTools::getRateDistribution(substory.getParams());
    }
    if (dynamic_cast<MixedSubstitutionModel*>(model))
      tl = new DRHomogeneousMixedTreeLikelihood(*tree, *sites, model, rDist, checkTree, true, true);
    else
      tl = new DRHomogeneousTreeLikelihood(*tree, *sites, model, rDist, checkTree);

    nbStates = model->getNumberOfStates();
  }
  else if (nhOpt == "one_per_branch")
  {
    model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, gCode.get(), sites, substory.getParams());
    if (model->getName() != "RE08") SiteContainerTools::changeGapsToUnknownCharacters(*sites);
    if (model->getNumberOfStates() > model->getAlphabet()->getSize())
    {
      //Markov-modulated Markov model!
      rDist = new ConstantRateDistribution();
    }
    else
    {
      rDist = PhylogeneticsApplicationTools::getRateDistribution(substory.getParams());
    }
    vector<double> rateFreqs;
    if (model->getNumberOfStates() != alphabet->getSize())
    {
      //Markov-Modulated Markov Model...
      unsigned int n =(unsigned int)(model->getNumberOfStates() / alphabet->getSize());
      rateFreqs = vector<double>(n, 1./(double)n); // Equal rates assumed for now, may be changed later (actually, in the most general case,
                                                   // we should assume a rate distribution for the root also!!!  
    }
    FrequenciesSet * rootFreqs = PhylogeneticsApplicationTools::getRootFrequenciesSet(alphabet, gCode.get(), sites, substory.getParams(), rateFreqs);
    vector<string> globalParameters = ApplicationTools::getVectorParameter<string>("nonhomogeneous_one_per_branch.shared_parameters", substory.getParams(), ',', "");
    modelSet = SubstitutionModelSetTools::createNonHomogeneousModelSet(model, rootFreqs, tree, globalParameters); 
    model = 0;
    if (dynamic_cast<MixedSubstitutionModelSet*>(modelSet))
      throw Exception("Non-homogeneous mixed substitution ancestor reconstruction not implemented, sorry!");
    tl = new DRNonHomogeneousTreeLikelihood(*tree, *sites, modelSet, rDist, true);
    nbStates = modelSet->getNumberOfStates();
  }
  else if (nhOpt == "general")
  {
    modelSet = PhylogeneticsApplicationTools::getSubstitutionModelSet(alphabet, gCode.get(), sites, substory.getParams());
    if (modelSet->getModel(0)->getName() != "RE08") SiteContainerTools::changeGapsToUnknownCharacters(*sites);
    if (modelSet->getNumberOfStates() > modelSet->getAlphabet()->getSize())
    {
      //Markov-modulated Markov model!
      rDist = new ConstantDistribution(1.);
    }
    else
    {
      rDist = PhylogeneticsApplicationTools::getRateDistribution(substory.getParams());
    }
    if (dynamic_cast<MixedSubstitutionModelSet*>(modelSet))
      throw Exception("Non-homogeneous mixed substitution ancestor reconstruction not implemented, sorry!");
    tl = new DRNonHomogeneousTreeLikelihood(*tree, *sites, modelSet, rDist, true);
    nbStates = modelSet->getNumberOfStates();
  }
  else throw Exception("Unknown option for nonhomogeneous: " + nhOpt);
  tl->initialize();
 
  delete tree;
    
  double logL = tl->getValue();
  if (isinf(logL))
  {
    // This may be due to null branch lengths, leading to null likelihood!
    ApplicationTools::displayWarning("!!! Warning!!! Likelihood is zero.");
    ApplicationTools::displayWarning("!!! This may be due to branch length == 0.");
    ApplicationTools::displayWarning("!!! All null branch lengths will be set to 0.000001.");
    ParameterList pl = tl->getBranchLengthsParameters();
    for(unsigned int i = 0; i < pl.size(); i++)
    {
      if(pl[i].getValue() < 0.000001) pl[i].setValue(0.000001);
    }
    tl->matchParametersValues(pl);
    logL = tl->getValue();
  }
  if (isinf(logL))
  {
    ApplicationTools::displayError("!!! Unexpected likelihood == 0.");
    ApplicationTools::displayError("!!! Looking at each site:");
    for(unsigned int i = 0; i < sites->getNumberOfSites(); i++)
    {
      (*ApplicationTools::error << "Site " << sites->getSite(i).getPosition() << "\tlog likelihood = " << tl->getLogLikelihoodForASite(i)).endLine();
    }
    ApplicationTools::displayError("!!! 0 values (inf in log) may be due to computer overflow, particularily if datasets are big (>~500 sequences).");
    exit(-1);
  }
  tree = new TreeTemplate<Node>(tl->getTree());

  // Write parameters to screen:
  ApplicationTools::displayResult("Log likelihood", TextTools::toString(tl->getValue(), 15));
  ParameterList parameters = tl->getSubstitutionModelParameters();
  for (unsigned int i = 0; i < parameters.size(); i++)
  {
    ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
  }
  parameters = tl->getRateDistributionParameters();
  for (unsigned int i = 0; i < parameters.size(); i++)
  {
    ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
  }

  // Getting posterior rate class distribution:
  DiscreteDistribution* prDist = RASTools::getPosteriorRateDistribution(*tl);
  ApplicationTools::displayMessage("\nPosterior rate distribution for dataset:\n");
  if (ApplicationTools::message) prDist->print(*ApplicationTools::message);
  ApplicationTools::displayMessage("\n");
  delete prDist;

  // Reconstruct ancestral sequences:
  string reconstruction = ApplicationTools::getStringParameter("asr.method", substory.getParams(), "marginal", "", true, false);
  ApplicationTools::displayResult("Ancestral state reconstruction method", reconstruction);
  bool probs = false;

  AncestralStateReconstruction *asr = 0;
  bool probMethod = false;
  if (reconstruction == "none")
  {
    //do nothing
  } else if (reconstruction == "marginal") {
    asr = new MarginalAncestralStateReconstruction(tl);
    probMethod = true;
  } else
    throw Exception("Unknown ancestral state reconstruction method: " + reconstruction);

  string outputFile;
  if (asr) {
    if (probMethod)
    {
      probs = ApplicationTools::getBooleanParameter("asr.probabilities", substory.getParams(), false, "", true, false);
      ApplicationTools::displayResult("Output probabilities", probs ? "yes" : "no");
    }

    // Write infos to file:
    outputFile = ApplicationTools::getAFilePath("output.sites.file", substory.getParams(), false, false);
    if (outputFile != "none")
    {
      ApplicationTools::displayResult("Output file for sites", outputFile);
      ofstream out(outputFile.c_str(), ios::out);
      TreeTemplate<Node> ttree(*tree);
      vector<Node *> nodes = ttree.getInnerNodes();
      size_t nbNodes = nodes.size();
    
      // Get the rate class with maximum posterior probability:
      vector<size_t> classes = tl->getRateClassWithMaxPostProbOfEachSite();
      // Get the posterior rate, i.e. rate averaged over all posterior probabilities:
      Vdouble rates = tl->getPosteriorRateOfEachSite();
      // Get the ancestral sequences:
      vector<Sequence*> sequences(nbNodes);
      vector<VVdouble*> probabilities(nbNodes);

      vector<string> colNames;
      colNames.push_back("Sites");
      colNames.push_back("is.complete");
      colNames.push_back("is.constant");
      colNames.push_back("lnL");
      colNames.push_back("rc");
      colNames.push_back("pr");
      for (size_t i = 0; i < nbNodes; i++) {
        Node *node = nodes[i];
        colNames.push_back("max." + TextTools::toString(node->getId()));
        if (probs) {
          probabilities[i] = new VVdouble();
          //The cast will have to be updated when more probabilistic method will be available:
          sequences[i] = dynamic_cast<MarginalAncestralStateReconstruction *>(asr)->getAncestralSequenceForNode(node->getId(), probabilities[i], false);

          for (unsigned int j = 0; j < nbStates; j++) {
            colNames.push_back("prob." + TextTools::toString(node->getId()) + "." + alphabet->intToChar((int)j));
          }
        }
        else
        {
          if (node->isLeaf()) {

          } else {
            sequences[i] = asr->getAncestralSequenceForNode(node->getId());
          }
        }
      }

      //Now fill the table:
      vector<string> row(colNames.size());
      DataTable* infos = new DataTable(colNames);
    
      for (size_t i = 0; i < sites->getNumberOfSites(); i++)
      {
        double lnL = tl->getLogLikelihoodForASite(i);
        const Site* currentSite = &sites->getSite(i);
        int currentSitePosition = currentSite->getPosition();
        string isCompl = "NA";
        string isConst = "NA";
        try { isCompl = (SiteTools::isComplete(*currentSite) ? "1" : "0"); }
        catch(EmptySiteException& ex) {}
        try { isConst = (SiteTools::isConstant(*currentSite) ? "1" : "0"); }
        catch(EmptySiteException& ex) {}
        row[0] = (string("[" + TextTools::toString(currentSitePosition) + "]"));
        row[1] = isCompl;
        row[2] = isConst;
        row[3] = TextTools::toString(lnL);
        row[4] = TextTools::toString(classes[i]);
        row[5] = TextTools::toString(rates[i]);

        unsigned int k = 6;
        for (unsigned int j = 0; j < nbNodes; j++) {
          row[k] = sequences[j]->getChar(i);
          k++;
          if (probs) {
            for (unsigned int l = 0; l < nbStates; l++) {
              row[k] = TextTools::toString((*probabilities[j])[i][l]);
              k++;
            }
          }
        }

        infos->addRow(row);
      }

      DataTable::write(*infos, out, "\t");

      delete infos;
    }

    SiteContainer* asSites = 0;
    if (probMethod)
    {
      bool sample = ApplicationTools::getBooleanParameter("asr.sample", substory.getParams(), false, "", true, false);
      ApplicationTools::displayResult("Sample from posterior distribution", sample ? "yes" : "no");
      if (sample)
      {
        unsigned int nbSamples = ApplicationTools::getParameter<unsigned int>("asr.sample.number", substory.getParams(), 1, "", true, false);
        asSites = new AlignedSequenceContainer(alphabet);
        for (unsigned int i = 0; i < nbSamples; i++)
        {
          ApplicationTools::displayGauge(i, nbSamples-1, '=');
          SequenceContainer *sampleSites = dynamic_cast<MarginalAncestralStateReconstruction *>(asr)->getAncestralSequences(true);
          vector<string> names = sampleSites->getSequencesNames();
          for (unsigned int j = 0; j < names.size(); j++)
            names[j] += "_" + TextTools::toString(i+1);
          sampleSites->setSequencesNames(names, false);
          SequenceContainerTools::append(*asSites, *sampleSites);
          delete sampleSites;
        }
        ApplicationTools::message->endLine();
      }
      else
      {
        asSites = asr->getAncestralSequences();
      }
    }
    else
    {
      asSites = asr->getAncestralSequences();
    }
  
    //Add existing sequence to output?
    bool addExtant = ApplicationTools::getBooleanParameter("asr.add_extant", substory.getParams(), false, "", true, false);
    if (addExtant) {
      SequenceContainerTools::append(*asSites, *sites);
    }

    //Write output:
    if (ApplicationTools::getStringParameter("output.sequence.file", substory.getParams(), "none") != "none") {
      SequenceApplicationTools::writeAlignmentFile(*asSites, substory.getParams());
    }
    delete asSites;

    delete asr;
  }

  outputFile = ApplicationTools::getAFilePath("output.nodes.file", substory.getParams(), false, false);
  if (outputFile != "none")
  {
    ApplicationTools::displayResult("Output file for nodes", outputFile);
    ofstream out(outputFile.c_str(), ios::out);

    //Add existing sequence to output?
    bool addExtant = ApplicationTools::getBooleanParameter("output.nodes.add_extant", substory.getParams(), false, "", true, false);
    
    map<int, vector<double> > frequencies;
    TreeLikelihoodTools::getAncestralFrequencies(*tl, frequencies, addExtant);
    
    vector<string> colNames;
    colNames.push_back("Nodes");
    for (unsigned int i = 0; i < tl->getNumberOfStates(); i++)
      colNames.push_back("exp" + tl->getAlphabetStateAsChar(i));
    for (unsigned int i = 0; i < tl->getNumberOfStates(); i++)
      colNames.push_back("eb" + tl->getAlphabetStateAsChar(i));

    //Now fill the table:
    vector<string> row(colNames.size());
    DataTable* infos = new DataTable(colNames);
    
    for (map<int, vector<double> >::iterator it = frequencies.begin(); it != frequencies.end(); it++)
    {
      row[0] = TextTools::toString(it->first);
      Vdouble ebFreqs = DRTreeLikelihoodTools::getPosteriorStateFrequencies(*tl, it->first);
      for (unsigned int i = 0; i < tl->getNumberOfStates(); i++)
      {
        row[i + 1] = TextTools::toString(it->second[i]);
      }
      for (unsigned int i = 0; i < tl->getNumberOfStates(); i++)
      {
        row[i + tl->getNumberOfStates() + 1] = TextTools::toString(ebFreqs[i]);
      }
      infos->addRow(row);
    }
    
    DataTable::write(*infos, out, "\t");

    delete infos;
  }

  //////////////////////////////////////////////////////////////////////////////
  //
  //   Ancestral sequence reconstruction over, now doing substitution mapping
  //
  //////////////////////////////////////////////////////////////////////////////
  
  
      //
    // Initialize the parameters for the mapping:
    //

    SubstitutionRegister* reg = 0;
    string regTypeDesc = ApplicationTools::getStringParameter("map.type", mapnh.getParams(), "All", "", true, false);
    string regType = "";
    map<string, string> regArgs;
    KeyvalTools::parseProcedure(regTypeDesc, regType, regArgs);
    bool stationarity = true;
    if (regType == "All")
    {
      stationarity = ApplicationTools::getBooleanParameter("stationarity", regArgs, true);
      reg = new ComprehensiveSubstitutionRegister(model, false);
    }
    else if (regType == "Total")
    {
      reg = new TotalSubstitutionRegister(model);
    }    
    else if (regType == "Selected"){  
      string subsList = ApplicationTools::getStringParameter("substitution.list", mapnh.getParams(), "All", "", true, false);
      reg = new SelectedSubstitutionRegister(model, subsList);  
    }
    else if (regType == "IntraAA")
    {
      if (AlphabetTools::isCodonAlphabet(alphabet))
      {
        reg = new AAInteriorSubstitutionRegister(dynamic_cast<CodonSubstitutionModel*>(model)); 
      }
      else
        throw Exception("Internal amino-acid categorization is only available for codon alphabet!");
    }
    else if (regType == "InterAA")
    {
      if (AlphabetTools::isCodonAlphabet(alphabet))
      {
        reg = new AAExteriorSubstitutionRegister(dynamic_cast<CodonSubstitutionModel*>(model)); 
      }
      else
        throw Exception("External amino-acid categorization is only available for codon alphabet!");
    }
    else if (regType == "GC")
    {
      stationarity = ApplicationTools::getBooleanParameter("stationarity", regArgs, true);
      if (AlphabetTools::isNucleicAlphabet(alphabet))
        reg = new GCSubstitutionRegister(dynamic_cast<NucleotideSubstitutionModel*>(model), false);
      else if (AlphabetTools::isCodonAlphabet(alphabet))
      {
        reg = new GCSynonymousSubstitutionRegister(dynamic_cast<CodonSubstitutionModel*>(model));
      }
      else
        throw Exception("GC categorization is only available for nucleotide or codon alphabets!");
    }
    else if (regType == "TsTv")
    {
      if (AlphabetTools::isNucleicAlphabet(alphabet))
        reg = new TsTvSubstitutionRegister(dynamic_cast<NucleotideSubstitutionModel*>(model));
      throw Exception("TsTv categorization is only available for nucleotide alphabet!");
    }

    else if (regType == "DnDs")
    {
      if (AlphabetTools::isCodonAlphabet(alphabet))
      {
        reg = new DnDsSubstitutionRegister(dynamic_cast<CodonSubstitutionModel*>(model), false);
      }
      else
        throw Exception("DnDs categorization is only available for codon alphabet!");
    }
    else
      throw Exception("Unsupported substitution categorization: " + regType);

    //Write categories:
    for (size_t i = 0; i < reg->getNumberOfSubstitutionTypes(); ++i)
      ApplicationTools::displayResult("  * Count type " + TextTools::toString(i + 1), reg->getTypeName(i + 1));

    // specific parameters to the null models
    string nullModelParams = ApplicationTools::getStringParameter("nullModelParams", mapnh.getParams(), "");

    ParameterList nullParams;
    if (nullModelParams != "")
    {
      string modelName = "";
      map<string, string> npv;
      KeyvalTools::multipleKeyvals(nullModelParams, npv, ",", false);

      map<string, string>::iterator mi(npv.begin());
      while (mi != npv.end())
      {
        nullParams.addParameter(Parameter(mi->first, TextTools::toDouble(mi->second)));
        ApplicationTools::displayResult("null Parameter " + mi->first, mi->second);
        
        mi++;
      }
    }

 
    //
    // Performs mapping
    //

    vector<int> ids = tl->getTree().getNodesId();
    ids.pop_back(); // remove root id.
    vector< vector<double> > counts;
    double thresholdSat = ApplicationTools::getDoubleParameter("count.max", mapnh.getParams(), -1);
    if (thresholdSat > 0)
      ApplicationTools::displayResult("Saturation threshold used", thresholdSat);

    if (nullModelParams != "")
    {
      if (model)
      {
        auto_ptr<SubstitutionModel> nullModel(model->clone());
        
        ParameterList pl;
        const ParameterList pl0 = nullModel->getParameters();
        
        for (size_t i = 0; i < nullParams.size(); ++i)
        {
          vector<string> pn = pl0.getMatchingParameterNames(nullParams[i].getName());
          for (size_t j = 0; j < pn.size(); ++j)
          {
            pl.addParameter(Parameter(pn[j], nullParams[i].getValue()));
          }
        }

        nullModel->matchParametersValues(pl);
        
        counts = SubstitutionMappingTools::getNormalizedCountsPerBranch(*tl, ids, model, nullModel.get(), *reg, true);
      }
      else
      {
        auto_ptr<SubstitutionModelSet> nullModelSet(modelSet->clone());
        ParameterList pl;
        const ParameterList pl0 = nullModelSet->getParameters();

        for (size_t i = 0; i < nullParams.size(); ++i)
        {
          vector<string> pn = pl0.getMatchingParameterNames(nullParams[i].getName());
          for (size_t j = 0; j < pn.size(); ++j)
          {
            pl.addParameter(Parameter(pn[j], nullParams[i].getValue()));
          }
        }

        nullModelSet->matchParametersValues(pl);

        counts = SubstitutionMappingTools::getNormalizedCountsPerBranch(*tl, ids, modelSet, nullModelSet.get(), *reg, true);
      }
    }
    else
      counts = SubstitutionMappingTools::getRelativeCountsPerBranch(*tl, ids, model ? model : modelSet->getModel(0), *reg, stationarity, thresholdSat);

    vector<string> outputDesc = ApplicationTools::getVectorParameter<string>("output.counts", mapnh.getParams(), ',', "PerType(prefix=)");
    for (vector<string>::iterator it = outputDesc.begin(); it != outputDesc.end(); ++it) {
      string outputType;
      map<string, string> outputArgs;
      KeyvalTools::parseProcedure(*it, outputType, outputArgs);
      if (outputType == "PerType")
      {
        // Write count trees:
        string treePathPrefix = ApplicationTools::getStringParameter("prefix", outputArgs, "mapping_counts_per_type_", "", true, 1);
        if (treePathPrefix != "none")
        {
          Newick newick;
          for (size_t i = 0; i < reg->getNumberOfSubstitutionTypes(); ++i)
          {
            string path = treePathPrefix + TextTools::toString(i + 1) + string(".dnd");
            ApplicationTools::displayResult(string("Output counts of type ") + TextTools::toString(i + 1) + string(" to file"), path);
            Tree* cTree = tree->clone();
            buildCountTree(counts, ids, cTree, i);
            newick.write(*cTree, path);
            delete cTree;
          }
        }
      }
      else if (outputType == "PerSite")
      {
        string perSitenf = ApplicationTools::getStringParameter("file", outputArgs, "mapping_counts_per_site.txt", "", true, 1);
        if (perSitenf != "none")
        {
          SubstitutionMappingTools::outputTotalCountsPerBranchPerSite(perSitenf, *tl, ids, model ? model : modelSet->getModel(0), *reg);
        }
      }
      else if (outputType == "PerSitePerType")
      {
        string tablePathPrefix = ApplicationTools::getStringParameter("prefix", outputArgs, "mapping_counts_per_site_per_type_", "", true, 1);
        if (tablePathPrefix != "none")
        {
          SubstitutionMappingTools::outputIndividualCountsPerBranchPerSite(tablePathPrefix, *tl, ids, model ? model : modelSet->getModel(0), *reg);
        }
      }
    }

    // Rounded counts
    vector< vector<size_t> > countsint;
    for (size_t i = 0; i < counts.size(); i++)
    {
      vector<size_t> countsi2;
      for (size_t j = 0; j < counts[i].size(); j++)
      {
        countsi2.push_back(static_cast<size_t>(floor( counts[i][j] + 0.5)));
      }
      countsint.push_back(countsi2);
    }


  

  delete alphabet;
  delete sites;
  if(model)    delete model;
  if(modelSet) delete modelSet;
  delete rDist;
  delete tl;
  delete tree;
  substory.done();

  }
  catch (exception & e)
  {
    cout << e.what() << endl;
    return 1;
  }

  return 0;
}




void buildCountTree(
  const vector< vector<double> >& counts,
  const vector<int>& ids,
  Tree* cTree,
  size_t type)
{
  for (size_t i = 0; i < ids.size(); ++i)
  {
    if (cTree->hasFather(ids[i]))
    {
      cTree->setDistanceToFather(ids[i], counts[i][type]);
    }
  }
}



/*
int main(int args, char** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*                     Map NH, version 0.2.0                      *" << endl;
  cout << "* Authors: J. Dutheil                       Created on  09/12/10 *" << endl;
  cout << "*          B. Boussau                       Modif. 17/12/11      *" << endl;
  cout << "*          L. Guéguen                       Last Modif. 17/06/13 *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  if (args == 1)
  {
    help();
    exit(0);
  }

  try
  {
    BppApplication mapnh(args, argv, "MapNH");
    mapnh.startTimer();

    Alphabet* alphabet = SequenceApplicationTools::getAlphabet(mapnh.getParams(), "", false);
    auto_ptr<GeneticCode> gCode;
    CodonAlphabet* codonAlphabet = dynamic_cast<CodonAlphabet*>(alphabet);
    if (codonAlphabet) {
      string codeDesc = ApplicationTools::getStringParameter("genetic_code", mapnh.getParams(), "Standard", "", true, true);
      ApplicationTools::displayResult("Genetic Code", codeDesc);
      
      gCode.reset(SequenceApplicationTools::getGeneticCode(codonAlphabet->getNucleicAlphabet(), codeDesc));
    }

    VectorSiteContainer* allSites = SequenceApplicationTools::getSiteContainer(alphabet, mapnh.getParams());
    VectorSiteContainer* sites = SequenceApplicationTools::getSitesToAnalyse(*allSites, mapnh.getParams());
    delete allSites;

    ApplicationTools::displayResult("Number of sequences", TextTools::toString(sites->getNumberOfSequences()));
    ApplicationTools::displayResult("Number of sites", TextTools::toString(sites->getNumberOfSites()));

    // Get the initial tree
    Tree* tree = PhylogeneticsApplicationTools::getTree(mapnh.getParams());
    ApplicationTools::displayResult("Number of leaves", TextTools::toString(tree->getNumberOfLeaves()));
    // Convert to NHX if input tree is newick or nexus?
    string treeIdOut = ApplicationTools::getAFilePath("output.tree_with_id.file", mapnh.getParams(), false, false);
    if (treeIdOut != "none")
    {
      Nhx nhx(true);
      nhx.write(*tree, treeIdOut);
    }

    //
    // Get substitution model and compute likelihood arrays
    //

    string nhOpt = ApplicationTools::getStringParameter("nonhomogeneous", mapnh.getParams(), "no", "", true, false);
    ApplicationTools::displayResult("Heterogeneous model", nhOpt);

    DRTreeLikelihood* tl     = 0;
    SubstitutionModel* model    = 0;
    SubstitutionModelSet* modelSet = 0;
    DiscreteDistribution* rDist    = 0;

    if (nhOpt == "no")
    {
      model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, gCode.get(), sites, mapnh.getParams());
      if (model->getName() != "RE08")
        SiteContainerTools::changeGapsToUnknownCharacters(*sites);
      if (model->getNumberOfStates() > model->getAlphabet()->getSize())
      {
        // Markov-modulated Markov model!
        rDist = new ConstantRateDistribution();
      }
      else
      {
        rDist = PhylogeneticsApplicationTools::getRateDistribution(mapnh.getParams());
      }
      
      tl = new DRHomogeneousTreeLikelihood(*tree, *sites, model, rDist, false, false);
    }
    else if (nhOpt == "one_per_branch")
    {
      model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, gCode.get(), sites, mapnh.getParams());
      if (model->getName() != "RE08")
        SiteContainerTools::changeGapsToUnknownCharacters(*sites);
      if (model->getNumberOfStates() > model->getAlphabet()->getSize())
      {
        // Markov-modulated Markov model!
        rDist = new ConstantRateDistribution();
      }
      else
      {
        rDist = PhylogeneticsApplicationTools::getRateDistribution(mapnh.getParams());
      }
      vector<double> rateFreqs;
      if (model->getNumberOfStates() != alphabet->getSize())
      {
        // Markov-Modulated Markov Model...
        size_t n = (size_t)(model->getNumberOfStates() / alphabet->getSize());
        rateFreqs = vector<double>(n, 1. / static_cast<double>(n)); // Equal rates assumed for now, may be changed later (actually, in the most general case,
        // we should assume a rate distribution for the root also!!!
      }
      FrequenciesSet* rootFreqs = PhylogeneticsApplicationTools::getRootFrequenciesSet(alphabet, gCode.get(), sites, mapnh.getParams(), rateFreqs);
      vector<string> globalParameters = ApplicationTools::getVectorParameter<string>("nonhomogeneous_one_per_branch.shared_parameters", mapnh.getParams(), ',', "");
      modelSet = SubstitutionModelSetTools::createNonHomogeneousModelSet(model, rootFreqs, tree, globalParameters);
      tl = new DRNonHomogeneousTreeLikelihood(*tree, *sites, modelSet, rDist, false, false);
    }
    else if (nhOpt == "general")
    {
      modelSet = PhylogeneticsApplicationTools::getSubstitutionModelSet(alphabet, gCode.get(), sites, mapnh.getParams());
      model = modelSet->getModel(0);
      if (modelSet->getModel(0)->getName() != "RE08")
        SiteContainerTools::changeGapsToUnknownCharacters(*sites);
      if (modelSet->getNumberOfStates() > modelSet->getAlphabet()->getSize())
      {
        // Markov-modulated Markov model!
        rDist = new ConstantDistribution(1.);
      }
      else
      {
        rDist = PhylogeneticsApplicationTools::getRateDistribution(mapnh.getParams());
      }
      tl = new DRNonHomogeneousTreeLikelihood(*tree, *sites, modelSet, rDist, false, false);
    }
    else
      throw Exception("Unknown option for nonhomogeneous: " + nhOpt);
    tl->initialize();

    ApplicationTools::displayResult("Log-Likelihood", tl->getLogLikelihood());

    //Check for saturation:
    double ll = tl->getValue();
    if (isinf(ll))
    {
      ApplicationTools::displayError("!!! Unexpected initial likelihood == 0.");
      if (codonAlphabet)
      {
        bool f = false;
        size_t s;
        for (size_t i = 0; i < sites->getNumberOfSites(); i++) {
          if (isinf(tl->getLogLikelihoodForASite(i))) {
            const Site& site = sites->getSite(i);
            s = site.size();
            for (size_t j = 0; j < s; j++) {
              if (gCode->isStop(site.getValue(j))) {
                (*ApplicationTools::error << "Stop Codon at site " << site.getPosition() << " in sequence " << sites->getSequence(j).getName()).endLine();
                f = true;
              }
            }
          }
        }
        if (f)
          exit(-1);
      }
      bool removeSaturated = ApplicationTools::getBooleanParameter("input.sequence.remove_saturated_sites", mapnh.getParams(), false, "", true, 1);
      if (!removeSaturated) {
        ofstream debug ("DEBUG_likelihoods.txt", ios::out);
        for (size_t i = 0; i < sites->getNumberOfSites(); i++)
        {
          debug << "Position " << sites->getSite(i).getPosition() << " = " << tl->getLogLikelihoodForASite(i) << endl; 
        }
        debug.close();
        ApplicationTools::displayError("!!! Site-specific likelihood have been written in file DEBUG_likelihoods.txt .");
        ApplicationTools::displayError("!!! 0 values (inf in log) may be due to computer overflow, particularily if datasets are big (>~500 sequences).");
        ApplicationTools::displayError("!!! You may want to try input.sequence.remove_saturated_sites = yes to ignore positions with likelihood 0.");
        exit(1);
      } else {
        for (size_t i = sites->getNumberOfSites(); i > 0; --i) {
          if (isinf(tl->getLogLikelihoodForASite(i - 1))) {
            ApplicationTools::displayResult("Ignore saturated site", sites->getSite(i - 1).getPosition());
            sites->deleteSite(i - 1);
          }
        }
        ApplicationTools::displayResult("Number of sites retained", sites->getNumberOfSites());
        tl->setData(*sites);
        tl->initialize();
        ll = tl->getValue();
        if (isinf(ll)) {
          throw Exception("Likelihood is still 0 after saturated sites are removed! Looks like a bug...");
        }
        ApplicationTools::displayResult("Initial log likelihood", TextTools::toString(-ll, 15));
      }
    }
    
 
    // Cleaning up:
    delete alphabet;
    delete sites;
    delete tree;
    if (modelSet)
      delete modelSet;
    else
      delete model;
    delete rDist;
    delete reg;
    mapnh.done();
  }
  catch (exception& e)
  {
    cout << e.what() << endl;
    exit(-1);
  }

  return 0;
}


*/
