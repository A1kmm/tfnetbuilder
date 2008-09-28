#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <iostream>
#include <fstream>
#include <boost/random.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

namespace po = boost::program_options;
namespace bll = boost::lambda;

int
main(int argc, char** argv)
{
  std::string model;

  po::options_description desc;

  desc.add_options()
    ("model", po::value<std::string>(&model), "TF net model to perturb")
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
  }

  if (wrong != "")
    std::cerr << "Missing option: " << wrong << std::endl;

  if (vm.count("help") || wrong != "")
  {
    std::cout << desc << std::endl;
    return 1;
  }

  std::ifstream modelFile(model.c_str());

  std::string l;
  std::getline(modelFile, l);
  if (l != "VERTICES")
  {
    std::cerr << "Expected VERTICES line" << std::endl;
    return 1;
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
      uint32_t n(strtoul(m[1].str().c_str(), NULL, 10));
      allNames.push_back(m[2]);
      uint32_t rv(rng());
      allNumbers.push_back(std::pair<uint32_t, uint32_t>(n, rv));
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

  return 0;
}
