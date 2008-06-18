#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include "../parsegenbank/GenbankParser.hpp"
#include <iostream>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/tokenizer.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/replace.hpp>

namespace po = boost::program_options;
namespace fs = boost::filesystem;
namespace io = boost::iostreams;

class Gene
{
public:
  Gene(uint32_t aOffset, uint32_t aHgncId = 1)
    : offset(aOffset), hgncId(aHgncId)
  {
  }

  bool
  operator<(const Gene& aGene) const
  {
    return (offset < aGene.offset);
  }

  uint32_t offset;
  uint32_t hgncId;
};

class TFNetBuilder
  : public GenBankSink
{
public:
  TFNetBuilder(const fs::path& aBaSeTraM)
    : mTFBSProcessed(0), mBaSeTraM(aBaSeTraM), mGBP(NewGenBankParser()),
      mBTP(NewGenBankParser()), mComplement(false), mTFBSSink(this)
  {
    mGBP->SetSink(this);
    mBTP->SetSink(&mTFBSSink);
  }

  ~TFNetBuilder()
  {
    delete mGBP;
    delete mBTP;
  }

  void
  generateOutput(std::ostream& aOutput)
  {
    aOutput << "VERTICES" << std::endl;
    for
    (
     std::set<uint32_t>::iterator i = mUsedHGNCIds.begin();
     i != mUsedHGNCIds.end();
     i++
    )
      aOutput << "VERTEX " << (*i) << " " << mNameByHGNCId[*i]
              << std::endl;
    aOutput << "ENDVERTICES" << std::endl;

    uint32_t nEdges = 0;

    // Now build a multimap of edges by source...
    std::multimap<uint32_t, uint32_t> edges;
    for
    (
     std::set<std::pair<uint32_t, uint32_t> >::iterator i = mEdges.begin();
     i != mEdges.end();
     i++
    )
    {
      nEdges++;
      edges.insert(*i);
    }

    for
    (
     std::set<uint32_t>::iterator i = mUsedHGNCIds.begin();
     i != mUsedHGNCIds.end();
     i++
    )
    {
      std::multimap<uint32_t, uint32_t>::iterator j(edges.find(*i));

      if (j == edges.end())
        continue;

      aOutput << "EDGES " << (*i)
              << " (";
      for
      (
       ;
       j != edges.end() && (*j).first == (*i);
       j++
      )
        aOutput << (*j).second << " ";
      aOutput << ")" << std::endl;
    }

    aOutput << "# There are " << nEdges << " edges" << std::endl;
    aOutput << "# " << mTFBSProcessed
            << " transcription factor binding sites processed."
            << std::endl;
  }

  void
  indexMatrices(const std::string& aPath)
  {
    static const boost::regex AcPat("^AC[ \\t]+(.*)$");
    static const boost::regex BfPat("^BF[ \\t]+[^ ]+ ([^;]*);.*$");
    static const boost::regex NaPat("^NA[ \\t]+([^ ]+).*$");

    io::filtering_istream db;
    db.push(io::file_source(aPath));

    bool seenAC(false);
    std::string AC;
    std::set<std::string> BF;

    while (db.good())
    {
      std::string line;
      std::getline(db, line);
      boost::smatch res;

      if (line == "//")
      {
        if (BF.size() && seenAC)
        {
          for
          (
           std::set<std::string>::iterator j(BF.begin());
           j != BF.end();
           j++
          )
          {
            // Look up the name from HGNC...
            std::map<std::string, uint32_t>::iterator i
              (mHGNCIdMappings.find(*j));
            

            if (i != mHGNCIdMappings.end())
              mHGNCByTRANSFAC.insert(std::pair<std::string, uint32_t>(AC, (*i).second));
          }
        }
        seenAC = false;
        BF.clear();
      }
      else if (boost::regex_match(line, res, AcPat))
      {
        AC = res[1];
        seenAC = true;
      }
      else if (boost::regex_match(line, res, BfPat))
      {
        BF.insert(cleanup_HGNC_name(res[1]));
      }
      else if (boost::regex_match(line, res, NaPat))
      {
        BF.insert(cleanup_HGNC_name(res[1]));
      }
    }
  }

  void
  loadHGNCDatabase(const std::string& aPath)
  {
    io::filtering_istream db;
    db.push(io::file_source(aPath));

    // Skip the header...
    std::string entry;
    std::getline(db, entry);

    while (db.good())
    {
      std::getline(db, entry);

      boost::tokenizer<boost::char_separator<char> >
        tok(entry, boost::char_separator<char>("\t", "",
                                               boost::keep_empty_tokens));
      std::vector<std::string> v(tok.begin(), tok.end());

      if (v.size() < 6)
        continue;

      if (v[3] != "Approved")
        continue;

      uint32_t hgncId = strtoul(v[0].c_str(), NULL, 10);
      addHGNCMapping(v[1], hgncId, true);

      static const boost::regex rtok("[, ]+");

      addHGNCMapping(v[2], hgncId, false);

      boost::sregex_token_iterator rti1
        (make_regex_token_iterator(v[4], rtok, -1));
      boost::sregex_token_iterator end;
      for (; rti1 != end; rti1++)
        addHGNCMapping(*rti1, hgncId, false);

      boost::sregex_token_iterator rti2
        (make_regex_token_iterator(v[5], rtok, -1));
      for (; rti2 != end; rti2++)
        addHGNCMapping(*rti2, hgncId, false);
    }
  }

  void
  processChromosome(const std::string& aFile)
  {
    TextSource* ts = NewBufferedFileSource(aFile.c_str());
    mGBP->SetSource(ts);

    mChromosomeDir = mBaSeTraM;
    mChromosomeDir /= fs::basename(aFile);

    try
    {
      mGBP->Parse();
    }
    catch (ParserException& pe)
    {
      std::cout << "Parse error: " << pe.what() << std::endl;
    }

    mGBP->SetSource(NULL);
    delete ts;
  }

  void
  OpenKeyword(const char* name, const char* value)
  {
    if (!::strcmp(name, "LOCUS"))
    {
      dealWithContig();
      
      mContigFile = mChromosomeDir;

      std::string locus(value);
      size_t pos = locus.find(" ");
      mContigFile /= locus.substr(0, pos);
    }
  }

  void
  CloseKeyword()
  {
  }

  void
  OpenFeature(const char* name, const char* location)
  {
    if (!::strcmp(name, "gene"))
    {
      if (!::strncmp(location, "complement(", 11))
      {
        mComplement = true;
        location += 11;
      }
      else
        mComplement = false;

      mGeneStart = strtoul(location, NULL, 10);

      location = strchr(location, '.');
      if (location == NULL)
        return;
      location += 2;

      mGeneEnd = strtoul(location, NULL, 10);
    }
  }

  void
  CloseFeature()
  {
  }

  void
  Qualifier(const char* name, const char* value)
  {
    if (strcmp(name, "db_xref"))
      return;

    if (strncmp(value, "HGNC:", 5))
      return;

    uint32_t hgncId = strtoul(value + 5, NULL, 10);

    // We now have a HGNC ID, a direction, and a start and end point.
    // Convert this to a range...

    if (mComplement)
      mReverseGenes.push_back(Gene(mGeneEnd, hgncId));
    else
      mForwardGenes.push_back(Gene(mGeneStart, hgncId));
  }

  void
  CodingData(const char* data)
  {
  }

  void
  dealWithContig()
  {
    if (mForwardGenes.size() == 0 && mReverseGenes.size() == 0)
      return;

    std::sort(mForwardGenes.begin(), mForwardGenes.end());
    std::sort(mReverseGenes.begin(), mReverseGenes.end());

    // Now we need to open the BaSeTraM output and start finding TFBSes...
    TextSource* ts = NewBufferedFileSource(mContigFile.string().c_str());
    mBTP->SetSource(ts);
    try
    {
      mBTP->Parse();
    }
    catch (const ParserException& pe)
    {
      std::cout << "Parse error: " << pe.what() << std::endl;
    }
    mBTP->SetSource(NULL);
    delete ts;

    mForwardGenes.clear();
    mReverseGenes.clear();
  }

  void
  processTFBS(bool isComplement, uint32_t start, uint32_t end,
              std::string TRANSFAC, double probability)
  {
    mTFBSProcessed++;
    if (isComplement)
    {
      size_t offset((start > kUpstreamZone) ? start - kUpstreamZone : 0);
      std::vector<Gene>::iterator next
        (std::upper_bound(mReverseGenes.begin(), mReverseGenes.end(),
                          Gene(start + kDownstreamZone)));

      for (next--;
           (next >= mReverseGenes.begin()) && (*next).offset >= offset;
           next--)
        processEdge(TRANSFAC, (*next).hgncId);
    }
    else
    {
      size_t offset((start > kDownstreamZone) ? start - kDownstreamZone : 0);
      std::vector<Gene>::iterator next
        (std::upper_bound(mForwardGenes.begin(), mForwardGenes.end(),
                          Gene(start + kUpstreamZone)));

      for (next--;
           (next >= mForwardGenes.begin()) && (*next).offset >= offset;
           next--)
        processEdge(TRANSFAC, (*next).hgncId);
    }
  }

private:
  uint32_t mTFBSProcessed;
  static const uint32_t kUpstreamZone = 200, kDownstreamZone = 50;
  fs::path mBaSeTraM, mChromosomeDir, mContigFile;
  GenBankParser* mGBP, * mBTP;
  bool mComplement;
  uint32_t mGeneStart, mGeneEnd;
  std::vector<Gene> mForwardGenes, mReverseGenes;

  std::map<std::string, uint32_t> mHGNCIdMappings, mHGNCByTRANSFAC;
  std::map<uint32_t, std::string> mNameByHGNCId;

  void addHGNCMapping(const std::string& aMapping, uint32_t aHGNC,
                      bool aOverride)
  {
    std::string dcmapping(cleanup_HGNC_name(aMapping));

    if (aOverride)
      mNameByHGNCId.insert(std::pair<uint32_t, std::string>(aHGNC, aMapping));

    std::map<std::string, uint32_t>::iterator i =
      mHGNCIdMappings.find(dcmapping);
    if (i != mHGNCIdMappings.end())
    {
      if (!aOverride)
        return;

      mHGNCIdMappings.erase(i);
    }

    mHGNCIdMappings.insert(std::pair<std::string, uint32_t>
                           (dcmapping, aHGNC));
  }

  void processEdge(const std::string& aTRANSFAC, uint32_t aTargetHGNC)
  {
    std::map<std::string, uint32_t>::iterator i
      (mHGNCByTRANSFAC.find(aTRANSFAC));
    if (i == mHGNCByTRANSFAC.end())
      return;

    uint32_t sourceHGNC = (*i).second;

    // We now have a source and target HGNC id... Just add them to the
    // edge set for now, and also mark the source and target as used.
    mUsedHGNCIds.insert(sourceHGNC);
    mUsedHGNCIds.insert(aTargetHGNC);
    mEdges.insert(std::pair<uint32_t, uint32_t>(sourceHGNC, aTargetHGNC));
  }

  std::string
  cleanup_HGNC_name(const std::string& aName)
  {
    std::string uc(boost::algorithm::to_upper_copy(aName));
    boost::algorithm::replace_all(uc, "-", "");

    return uc;
  }

  std::set<uint32_t> mUsedHGNCIds;
  std::set<std::pair<uint32_t, uint32_t> > mEdges;

  class TFBSSink
    : public GenBankSink
  {
  public:
    TFBSSink(TFNetBuilder* aTFBuilder)
      : mTFBuilder(aTFBuilder)
    {
    }

    void
    OpenKeyword(const char* name, const char* value)
    {
    }
    
    void
    CloseKeyword()
    {
    }
    
    void
    OpenFeature(const char* name, const char* location)
    {
      if (strcmp(name, "TFBS"))
      {
        mInTFBS = false;
        return;
      }

      mInTFBS = true;

      if (!strncmp(location, "complement(", 11))
      {
        mIsComplement = true;
        location += 11;
      }
      else
        mIsComplement = false;

      char* p;
      mStart = strtoul(location, &p, 10);
      p += 2;
      mEnd = strtoul(p, NULL, 10);
    }

    void
    CloseFeature()
    {
      if (!mInTFBS)
        return;

      mInTFBS = false;

      mTFBuilder->processTFBS(mIsComplement, mStart, mEnd, mTRANSFAC,
                              mProbability);
    }
    
    void
    Qualifier(const char* name, const char* value)
    {
      if (!strcmp(name, "probability"))
        mProbability = strtod(value, NULL);
      else if (!strcmp(name, "db_xref") &&
               !strncmp(value, "TRANSFAC:", 9))
        mTRANSFAC = value + 9;
    }

    void
    CodingData(const char* data)
    {
    }
  private:
    TFNetBuilder* mTFBuilder;
    std::string mTRANSFAC;
    double mProbability;
    bool mInTFBS, mIsComplement;
    uint32_t mStart, mEnd;
  };

  TFBSSink mTFBSSink;
};

int
main(int argc, char** argv)
{
  std::string basetram, genbank, hgnc, matrices;

  po::options_description desc;

  desc.add_options()
    ("basetram", po::value<std::string>(&basetram), "Location of BaSeTraM output directory")
    ("genbank", po::value<std::string>(&genbank), "Directory containing GenBank files")
    ("hgnc", po::value<std::string>(&hgnc), "File containing the HGNC names database")
    ("matrices", po::value<std::string>(&matrices), "File containing the TRANSFAC matrices "
     "database")
    ("help", "produce help message")
    ;
  
  po::variables_map vm;

  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  std::string wrong;
  if (!vm.count("help"))
  {
    if (!vm.count("basetram"))
      wrong = "basetram";
    else if (!vm.count("genbank"))
      wrong = "genbank";
    else if (!vm.count("hgnc"))
      wrong = "hgnc";
    else if (!vm.count("matrices"))
      wrong = "matrices";
  }

  if (wrong != "")
    std::cerr << "Missing option: " << wrong << std::endl;
  if (vm.count("help") || wrong != "")
  {
    std::cout << desc << std::endl;
    return 1;
  }

  if (!fs::is_directory(basetram))
  {
    std::cerr << "Supplied BaSeTraM 'directory' is not a valid directory."
              << std::endl;
    return 1;
  }

  if (!fs::is_directory(genbank))
  {
    std::cerr << "Supplied GenBank 'directory' is not a valid directory."
              << std::endl;
    return 1;
  }

  TFNetBuilder tfnb(basetram);

  // Now we start iterating through the GenBank files...
  for (fs::directory_iterator it(genbank); it != fs::directory_iterator(); it++)
  {
    if (fs::extension(it->path()) != ".gbk")
      continue;

    tfnb.loadHGNCDatabase(hgnc);
    tfnb.indexMatrices(matrices);
    try
    {
      tfnb.processChromosome(it->path().string());
    }
    catch (const ParserException& pe)
    {
      std::cout << "Parser error: " << pe.what() << std::endl;
    }
    tfnb.dealWithContig();
  }

  tfnb.generateOutput(std::cout);
}
