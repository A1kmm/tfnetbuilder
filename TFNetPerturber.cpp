#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <iostream>
#include <fstream>
#include <boost/random.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <map>

namespace po = boost::program_options;
namespace bll = boost::lambda;

class ModelPerturber;

class ModelPerturber
{
public:
  ModelPerturber(const char* aName)
    : mName(aName)
  {
    sRegistry.insert(std::pair<std::string, ModelPerturber*>(aName, this));
  }

  static ModelPerturber* findPerturberByName(const std::string& aName)
  {
    RegistryType::iterator i = sRegistry.find(aName);
    if (i == sRegistry.end())
      return NULL;
    
    return (*i).second;
  }

  static void listAvailablePerturbers(std::ostream& aOut)
  {
    aOut << "Available perturbers:" << std::endl;
    for (RegistryType::iterator i = sRegistry.begin();
         i != sRegistry.end();
         i++)
      aOut << "\t* " << (*i).first << std::endl
           << "\t\t* Parameter choices: " << (*i).second->getParameterHelp()
           << std::endl;
  }

  virtual const char* getParameterHelp() = 0;
  virtual void perturb(const std::string& aModelFile) = 0;
  const char* name()
  {
    return mName;
  }

  virtual void setParams(const std::string& aParams)
  {
  }

private:
  const char* mName;
  typedef std::map<std::string, ModelPerturber*> RegistryType;
  static RegistryType sRegistry;
};

ModelPerturber::RegistryType ModelPerturber::sRegistry;

class LabelSwitchingPerturber
  : public ModelPerturber
{
public:
  LabelSwitchingPerturber()
    : ModelPerturber("label_switching"), mProb(1.0)
  {
  }

  void
  setParams(const std::string& aParams)
  {
    mProb = strtod(aParams.c_str(), NULL);
  }

  void
  perturb(const std::string& aModelFile)
  {
    std::ifstream modelFile(aModelFile.c_str());

    std::string l;
    std::getline(modelFile, l);
    if (l != "VERTICES")
    {
      std::cerr << "Expected VERTICES line" << std::endl;
      return;
    }
    std::cout << l << std::endl;
    
    static const boost::regex vertexr("^VERTEX ([0-9]+) (.*)$");
    
    std::list<std::string> allNames;
    std::vector<std::pair<uint32_t, uint32_t> > allNumbers;
    boost::mt19937 rng;
    
    rng.seed(time(0));
    
    while (modelFile.good())
    {
      std::getline(modelFile, l);
      if (l == "ENDVERTICES")
        break;
      
      boost::smatch m;
      if (boost::regex_match(l, m, vertexr))
      {
        boost::uniform_real<double> ur;
        // Include it in the shuffle pool with probability mProb.
        if (ur(rng) < mProb)
        {
          uint32_t n(strtoul(m[1].str().c_str(), NULL, 10));
          allNames.push_back(m[2]);
          uint32_t rv(rng());
          allNumbers.push_back(std::pair<uint32_t, uint32_t>(n, rv));
        }
        else
        {
          std::cout << "VERTEX " << m[1].str() << m[2].str() << std::endl;
        }
      }
    }

    // Now scramble allNumbers...
    std::sort(
              allNumbers.begin(), allNumbers.end(),
              bll::bind<uint32_t>(&std::pair<uint32_t, uint32_t>::second, bll::_1) <
              bll::bind<uint32_t>(&std::pair<uint32_t, uint32_t>::second, bll::_2)
             );

    // Write out the vertices with their new names...
    std::list<std::string>::iterator i;
    std::vector<std::pair<uint32_t, uint32_t> >::iterator j;
    for (i = allNames.begin(), j = allNumbers.begin(); i != allNames.end();
         i++, j++)
    {
      std::cout << "VERTEX " << (*j).first << " " << (*i) << std::endl;
    }
    std::cout << "ENDVERTICES" << std::endl;

    while (modelFile.good())
    {
      std::getline(modelFile, l);
      std::cout << l << std::endl;
    }
  }

  const char* getParameterHelp()
  {
    return "Use --params=<probability> to specify the probability a genes / label is "
      "included in the pool to be scrambled.";
  }
  
private:
  double mProb;
};

class EdgeDeletingPerturber
  : public ModelPerturber
{
public:
  EdgeDeletingPerturber()
    : ModelPerturber("edge_deleting")
  {
  }

  void
  setParams(const std::string& aParams)
  {
  }

  void
  perturb(const std::string& aModelFile)
  {
    std::ifstream modelFile(aModelFile.c_str());

    std::string l;
    std::getline(modelFile, l);
    if (l != "VERTICES")
    {
      std::cerr << "Expected VERTICES line" << std::endl;
      return;
    }
    std::cout << l << std::endl;
  }
};


int
main(int argc, char** argv)
{
  std::string model, type, params;

  po::options_description desc;

  desc.add_options()
    ("model", po::value<std::string>(&model), "TF net model to perturb")
    ("type", po::value<std::string>(&type), "Type of perturber to use. --type=help to list")
    ("params", po::value<std::string>(&params), "Parameters for the perturber (type dependent)")
    ("help", "produce help message")
    ;

  po::variables_map vm;

  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  std::string wrong;
  if (!vm.count("help"))
  {
    if (!vm.count("model"))
      wrong = "model";
    else if (!vm.count("type"))
      wrong = "type";
  }

  if (wrong != "")
    std::cerr << "Missing option: " << wrong << std::endl;

  if (vm.count("help") || wrong != "")
  {
    std::cerr << desc << std::endl;
    return 1;
  }

  if (type == "help")
  {
    ModelPerturber::listAvailablePerturbers(std::cout);
    return 1;
  }

  ModelPerturber* mp = ModelPerturber::findPerturberByName(type);
  if (mp == NULL)
  {
    std::cerr << "Invalid model perturber type requested."
              << std::endl;
    return 1;
  }

  if (vm.count("params"))
    mp->setParams(params);
  mp->perturb(model);

  return 0;
}
