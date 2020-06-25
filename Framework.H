#ifndef _FRAMEWORK_
#define _FRAMEWORK_
#include <map>
#include <vector>
#include <string>
#include <stdio.h>
using namespace std;
class GFFReader;
class GeneNameMapper;
class Framework
{
	public:
	Framework();
	~Framework();
	int readGFFFile(const char*);
	int readChromMap(const char*);
	int readSites(const char*,const char*);
        int mapSitesToGenes(const char*,const char*,int windowSizeUp=2000,int windowSizeDown=2000);
	int readRefFasta(const char* aFName,const char*);
	
	private:
	int readSitesSingle(const char*);
	int readSitesDouble(const char*);
	struct SiteInfo
	{
		int pos;
		int length;
		string chrom;
		char snpClass;
		string name;
	};
	vector<SiteInfo*> siteSet;
	GFFReader gffreader;
	int getChromosomeLength(string&);

	map<string,map<int,int>*> geneSiteMap;
	map<string,int> refchromLength;
	GeneNameMapper gnm;
	map<string,string> srcTargetChromMap;
};
#endif
