#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include "Gene.H"
#include "GFFReader.H"
#include "GeneNameMapper.H"
#include "Framework.H"

Framework::Framework()
{
	gnm.readGeneNames();
}

Framework::~Framework()
{
}


int 
Framework::readChromMap(const char* aFName)
{	
	ifstream inFile(aFName);
	char buffer[1024];
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		string srcChrom;
		string targetChrom;
		char*  tok=strtok(buffer,"\t");
		int tokCnt=0;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				srcChrom.append(tok);
			}
			else if(tokCnt==1)
			{
				targetChrom.append(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		srcTargetChromMap[srcChrom]=targetChrom;
	}
	inFile.close();
	return 0;
}

int
Framework::readSites(const char* aFName, const char* type)
{
	if(strcmp(type,"single")==0)
	{
		readSitesSingle(aFName);
	}
	else if(strcmp(type,"double")==0)
	{
		readSitesDouble(aFName);
	}
	return 0;
}

int 
Framework::readSitesSingle(const char* aFName)
{
	ifstream inFile(aFName);
	char buffer[1024];
	int lineCnt=0;
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		string chr;
		int pos=-1;
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		Framework::SiteInfo*  site=new Framework::SiteInfo;
		while(tok!=NULL)
		{
			switch(tokCnt)
			{
				case 0:
				{
					string srcKey(tok);
					if(srcTargetChromMap.find(srcKey)==srcTargetChromMap.end())
					{
						cout <<"No chromosome found for " << srcKey<< endl;
					}	
					site->chrom.append(srcTargetChromMap[srcKey]);
					break;
				}
				case 1:
				{
					site->pos=atoi(tok);
					break;
				}
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		siteSet.push_back(site);
	}
	inFile.close();
	return 0;
}


int 
Framework::readSitesDouble(const char* aFName)
{
	ifstream inFile(aFName);
	char buffer[1024];
	int lineCnt=0;
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		string chr;
		int pos=-1;
		int begin=0;
		int end=0;
		string name;
		char* tok=strtok(buffer," \t");
		int tokCnt=0;
		Framework::SiteInfo*  site=new Framework::SiteInfo;
		while(tok!=NULL)
		{
			char* tempTok=tok;
			char* pos=strchr(tok,'"');
			if(pos!=NULL)
			{
				tempTok=tok+1;
			}
			char* endpos=strchr(tempTok,'"');
			if(endpos!=NULL)
			{
				*endpos='\0';
			}
			switch(tokCnt)
			{
				case 0:
				{
					string srcKey(tempTok);
					if(srcTargetChromMap.find(srcKey)==srcTargetChromMap.end())
					{
						cout <<"No chromosome found for " << srcKey<< endl;
					}	
					else
					{
						site->chrom.append(srcTargetChromMap[srcKey]);
					}
					break;
				}
				case 1:
				{
					begin=atoi(tempTok);
					break;
				}
				case 2:
				{
					end=atoi(tempTok);
					break;
				}
				case 3:
				{
					name.append(tok);
				}		
			}
			tok=strtok(NULL," \t");
			tokCnt++;
		}
		if(site->chrom.length()==0)
		{
			delete site;
			continue;
		}
		pos=begin+(end-begin)/2;
		site->pos=pos;
		int length=(end-begin);
		site->length=length;
		site->name=name;
		siteSet.push_back(site);
	}
	inFile.close();
	return 0;
}

int 
Framework::readGFFFile(const char* aFName)
{
	gffreader.readGFFFile(aFName);
	gffreader.setNeighbors();
	return 0;
}

int
Framework::readRefFasta(const char* aFName, const char* chrnamemapping)
{
	char buffer[1024];
	ifstream cmapFile(chrnamemapping);
	map<string,string> nameIDMap;
	while(cmapFile.good())
	{
		cmapFile.getline(buffer,1023);
		char* tok=strtok(buffer,"\t");
		int tokcnt=0;
		string key1;
		string key2;
		while(tok!=NULL)
		{
			if(tokcnt==0)
			{
				key1.append(tok);		
			}
			else if(tokcnt==1)
			{
				key2.append(tok);
			}
			tok=strtok(NULL,"\t");
			tokcnt++;
		}
		nameIDMap[key1]=key2;
	}
	string bases;
	string chromosome;
	ifstream inFile(aFName);
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		if(strstr(buffer,">")!=NULL)
		{
			if(bases.length()>0)
			{
				refchromLength[chromosome]=bases.length();
			}
			bases.clear();
			chromosome.clear();
			char* tok=strchr(buffer,' ');
			string key;
			if(tok!=NULL)
			{
				*tok='\0';
			}
			key.append(buffer+1);
			if(nameIDMap.find(key)==nameIDMap.end())
			{
				cout <<"No chromosome " << key << endl;
				exit(0);
			}
			chromosome.append(nameIDMap[key]);
		}
		else
		{
			bases.append(buffer);
		}
	}
	if(bases.length()>0)
	{
		refchromLength[chromosome]=bases.length();
	}
	return 0;
}

int
Framework::getChromosomeLength(string& key)
{
	if(refchromLength.find(key)==refchromLength.end())
	{
		cout <<"No chromosome by name " << key << endl << endl << endl;
		//exit(0);
	}
	int len=refchromLength[key];
	return len;
}

int 
Framework::mapSitesToGenes(const char* outFName,const char* pFName,int windowSizeUp,int windowSizeDown)
{
	ofstream oFile(outFName);
	ofstream pFile(pFName);
	//Will take this as a parameter too.
	//int windowSize=2000;
	for(int s=0;s<siteSet.size();s++)
	{
		Framework::SiteInfo* site=siteSet[s];
		map<string,Gene*>* genesOnChrom=gffreader.getGeneSetForChromosome(site->chrom);
		//cout << site->chrom << "\t" << site->pos << "\t" << genesOnChrom->size() << endl;
		if(genesOnChrom==NULL)
		{
		        cout << "No chromosome " << site->chrom << " in site "<< site->pos << endl << endl << endl ;
			exit(0);
		}
		Gene* closestGene=NULL;
		int minDist=100000;
		char possibleHitType='\0';
		int position=site->pos;
		int hits=0;
		for(map<string,Gene*>::iterator gIter=genesOnChrom->begin();gIter!=genesOnChrom->end();gIter++)
		{
			Gene* gene=gIter->second;
			//cout << "Testing " << gene->name.c_str() << " at " << gene->geneStart << endl;
			int start=gIter->second->geneStart;
			int end=gIter->second->geneEnd;
			bool hitReg=false;
			int dist=0;
			if(gene->strand=='-')
			{
				//Regulatory region is +n bp or -n bp of the ATG
				int beginRegRegion=start-windowSizeDown;
				int endRegRegion=start+windowSizeUp;
				if(position>=beginRegRegion && position<=endRegRegion)
				{
					hitReg=true;
					hits++;
					dist=position-start;
				}
			}
			else 
			{
				int beginRegRegion=start-windowSizeUp;
				int endRegRegion=start+windowSizeDown;
				if(position>=beginRegRegion && position<=endRegRegion)
				{
					hitReg=true;
					hits++;
					dist=start-position;
				}
			}
			if(hitReg )
			{
				cout <<"Mapping " << site->chrom <<":" << site->pos << " to " << gene->name.c_str() << "\t" << start << "\t" << dist << endl;
				oFile << site->name.c_str() << "\t" << gene->name.c_str() << "\t" << site->chrom << "\t" << site->pos-(site->length/2) << "\t" << site->pos+(site->length/2) << "\t"<< dist << endl;			
			}
		}
		if(hits==0)
		{
			if(closestGene!=NULL)
			{
				cout <<"No hits found for site " << s <<  " at " << site->pos << " closest_gene " << closestGene->name << "\tDist="  << minDist
				<< closestGene->strand << "\t" << closestGene->geneStart << "\t" << closestGene->geneEnd<<  "\t" << possibleHitType << endl;
			}
			else
			{
				cout <<"No hits found for site " <<  s<< " closest gene null " << endl;
			}
			pFile << site->chrom << "\t" << site->pos-(site->length/2) << "\t" << site->pos+(site->length/2) << endl;
		}
	}
	/*ofstream oFile(outFName);
	for(map<string,map<int,int>*>::iterator gIter=geneSiteMap.begin();gIter!=geneSiteMap.end();gIter++)
	{
		//oFile << gnm.getCommonName(gIter->first.c_str()) << "\t" << gIter->first<<"\t" << gIter->second->size() << endl;
	  oFile << (string)gIter->first.c_str() << "\t"  << gIter->second->size() << endl;
	}*/
	oFile.close();
	pFile.close();
}

int
main(int argc, const char** argv)
{
	if(argc!=9)
	{
		cout <<"Usage: ./mapPeaksToGenes gfffile chromnamemap sites single|double output windowSizeUp windowSizeDown pFName"<< endl;
		return 0;
	}
	Framework fw;
	fw.readGFFFile(argv[1]);
	fw.readChromMap(argv[2]);
	fw.readSites(argv[3],argv[4]);
	int WindowSizeUp = atoi(argv[6]);
	int WindowSizeDown = atoi(argv[7]);	
	fw.mapSitesToGenes(argv[5],argv[8],WindowSizeUp,WindowSizeDown);
	return 0;
}
