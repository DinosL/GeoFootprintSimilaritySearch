#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include "bst.h"
#include "RTree.h"
#include <math.h>
#include <time.h>
#include <iostream>

#define LEN(arr) ((int) (sizeof (arr) / sizeof (arr)[0]))

#define SIMILARITY 0
#define SIMILARITY_AND_NORM 1
#define SPATIAL_SIMILARITY 2

using namespace std;

#define DIMENSIONS 2

typedef struct Triple {
	float value;
	int id;
	int type;  // 1 = Start, 0 = end
} Triple;


typedef struct Quad {
	float value;
	int id;
	int source; // 1 = S, 0 = R
	int type;  // 1 = Start, 0 = end
} Quad;

struct Rect
{
  Rect()  {}

  Rect(float a_minX, float a_minY, float a_maxX, float a_maxY)
  {
    min[0] = a_minX;
    min[1] = a_minY;

    max[0] = a_maxX;
    max[1] = a_maxY;
  }


  float min[2];
  float max[2];
};



int compare(const void * a, const void * b)
{

  Triple *tripleA = (Triple *)a;
  Triple *tripleB = (Triple *)b;

  if(tripleA->value > tripleB->value)
  	return 1;
  else if (tripleA->value < tripleB->value)
  	return -1;
  else
  	if(tripleA->type > tripleB->type)
  		return -1;
  	else
  		return 1;
}

int compareQuads(const void * a, const void * b)
{

  Quad *quadA = (Quad *)a;
  Quad *quadB = (Quad *)b;

  if(quadA->value > quadB->value)
  	return 1;
  else if (quadA->value < quadB->value)
  	return -1;
  else
  	if(quadA->type >= quadB->type)
  		return -1;
  	else
  		return 1;
}

int cmp (const void * a, const void * b)
{
  float fa = *(const float*) a;
  float fb = *(const float*) b;
  return (fa > fb) - (fa < fb);
}


int readData(FILE *f, int **userRegionsSize, Region ***userRegions, int *numUsers)
{
	char *line = NULL; // used for fileread
	size_t len = 0; // used for fileread
	ssize_t read; // used for fileread
	const char delim[2] = " "; // used for fileread
	int curuser = -1;
	int curday = -1;
	int curcount = 0;
	int curst; int curend;
	float curmbr[2*DIMENSIONS];
	int user, day;
	float x1,y1,x2,y2;
	char *token;
	int *countUserRegions;
	int i,j, currLength=-1;


	read = getline(&line, &len, f);
	if (read==-1)
	{
		printf("ERROR: first line is empty. Exiting...\n");
		return -1;
	};	
	token = strtok(line,delim);
	*numUsers = atoi(token);

	countUserRegions =  (int *)calloc(*numUsers,sizeof(int));	
	*userRegionsSize = (int *) calloc(*numUsers,sizeof(int));	

	// printf("numUsers = %d\n", *numUsers);
	int id = -1;
	int prevUser = -1;
	while ((read = getline(&line,&len,f)) != -1)	{

		token = strtok(line,delim);
		user = atoi(token);
		if(user != prevUser)
			id++;
		//Skip 4 items
		token = strtok(NULL,delim);
		token = strtok(NULL,delim);
		token = strtok(NULL,delim);
		token = strtok(NULL,delim);

		token = strtok(NULL,delim);
		x1 = atof(token);

		token = strtok(NULL,delim);
		y1 = atof(token);

		token = strtok(NULL,delim);
		x2 = atof(token);

		token = strtok(NULL,delim);
		y2 = atof(token);

		countUserRegions[id]+=1;
		prevUser = user;
	}
	*numUsers = id+1;				// It is id+1 instead of id. That way for loops dont need = to go up to actual user id

	*userRegions = (Region **) malloc(*numUsers*sizeof(Region *));
	for(i=0;i<*numUsers;i++)
	{
		(*userRegions)[i] = (Region *) malloc(countUserRegions[i]*sizeof(Region));
	}

	fseek(f, 0, SEEK_SET);
	read = getline(&line, &len, f);
	id = -1;
	prevUser = -1;
	while ((read = getline(&line,&len,f)) != -1)	{

		token = strtok(line,delim);
		user = atoi(token);
		if(user != prevUser)
			id++;
		// if(user>0)
		// 	break;
		//Skip 4 items
		token = strtok(NULL,delim);
		token = strtok(NULL,delim);
		token = strtok(NULL,delim);
		token = strtok(NULL,delim);

		token = strtok(NULL,delim);
		x1 = atof(token);

		token = strtok(NULL,delim);
		y1 = atof(token);

		token = strtok(NULL,delim);
		x2 = atof(token);

		token = strtok(NULL,delim);
		y2 = atof(token);

		currLength = (*userRegionsSize)[id];
		(*userRegions)[id][currLength].xStart = x1;
		(*userRegions)[id][currLength].yStart = y1;
        (*userRegions)[id][currLength].xEnd = x2;
        (*userRegions)[id][currLength].yEnd = y2;
        (*userRegionsSize)[id]++; 
		prevUser = user;
	}

    free(countUserRegions);
	return 0;
}

void computeNsqInorder(struct node* node, float *nsq, float prev, float tripleValue, struct node* root)
{
    if (node == NULL)
        return;
 
    /* first recur on left child */
    computeNsqInorder(node->left, nsq, prev, tripleValue,root);
    struct node* next = findSuccessor(root,node);
    if(next != NULL)
    {
    	*nsq += (next->value - node->value)*(tripleValue-prev)*node->frequency*node->frequency;
    }

    /* now recur on right child */
    computeNsqInorder(node->right, nsq, prev, tripleValue,root);
}


struct node *updateDStructure(struct node *D, struct Triple triple, Region *userRegions)
{
	if(triple.type == 1)			// START
	{
		struct node *e = findMaxforN(D,userRegions[triple.id].yStart);
		D = insert(D,userRegions[triple.id].yStart,e->frequency+1);
		if (e->value < userRegions[triple.id].yStart)
        {
            e = findSuccessor(D,e);
        }
        
		struct node *next = findSuccessor(D,e);
		
		while(next != NULL && next->value < userRegions[triple.id].yEnd)
		{
			next->frequency += 1;
			e = next;
			next = findSuccessor(D,e);
		}
		
        D = insert(D,userRegions[triple.id].yEnd,e->frequency-1);
	}
	else
	{
		struct node *e = find(D,userRegions[triple.id].yStart);

		struct node *next = findSuccessor(D,e);
		float v1 = e->value;

		while(next != NULL && next->value < userRegions[triple.id].yEnd)
		{
			if(next->frequency > 0)
				next->frequency -= 1;
			e = next;
			next = findSuccessor(D,e);
		}



        if (v1 == userRegions[triple.id].yEnd)
        {
            D = deleteNode(D,v1);
            
        }
        else if(next!=NULL)
		{
			D = deleteNode(D,next->value);
		}

         D = deleteNode(D,v1);
        
	}
	return D;
}


float computeNorm(int userRegionsSize, Region *userRegions)
{
	int j,i,k;
	float nsq, prev;


	struct Triple triples[2* userRegionsSize];
	k=0;
	for(j=0;j<userRegionsSize;j+=1)   // Use x-dimension
	{

			struct Triple startTemp = {userRegions[j].xStart,j,1};  // Add start of interval
			triples[k] = startTemp;
			k++;
			struct Triple endTemp = {userRegions[j].xEnd,j,0};  // Add end of interval
			triples[k] = endTemp;
			k++;
		
	}
	

	qsort(triples,k,sizeof(Triple),compare);
	
	nsq = 0.0;
	prev = triples[0].value;
	struct node *D = NULL;
	D = insert(D, 0, 0);
	
	for(j=0;j<k;j++)
	{
		computeNsqInorder(D, &nsq, prev, triples[j].value,D);
		D = updateDStructure(D, triples[j], userRegions);
		prev = triples[j].value;
	}
	return sqrt(nsq);
}


struct node* findMin(struct node* a, struct node* b)
{
	if(a->value > b->value)
		return b;
	else
		return a;
}

struct node *updateDStructureQuads(struct node *D, struct Quad quad, Region *userRegions)
{
	if(quad.type == 1)			// START
	{
		struct node *e = findMaxforN(D,userRegions[quad.id].yStart);
		D = insert(D,userRegions[quad.id].yStart,e->frequency+1);
		if (e->value < userRegions[quad.id].yStart)
        {
            e = findSuccessor(D,e);
        }
        
		struct node *next = findSuccessor(D,e);
		
		while(next != NULL && next->value < userRegions[quad.id].yEnd)
		{
			next->frequency += 1;
			e = next;
			next = findSuccessor(D,e);
		}
		
        D = insert(D,userRegions[quad.id].yEnd,e->frequency-1);
	}
	else
	{
		struct node *e = find(D,userRegions[quad.id].yStart);
		struct node *next = findSuccessor(D,e);
		float v1 = e->value;

		while(next != NULL && next->value < userRegions[quad.id].yEnd)
		{
			if(next->frequency > 0)
				next->frequency -= 1;
			e = next;
			next = findSuccessor(D,e);
		}

        if (v1 == userRegions[quad.id].yEnd)
        {
            D = deleteNode(D,v1);
            
        }
        else if(next!=NULL)
		{
			D = deleteNode(D,next->value);
		}

         D = deleteNode(D,v1);
        
	}
	return D;
}


float computeSimilarity(int userRegionsSizeR, Region *userRegionsR, int userRegionsSizeS, Region *userRegionsS, float normrsq, float normssq, int user1, int user2)
{
	int i,j,k;
	float prev, simn;
	struct node *up;

	int total_size = 2*(userRegionsSizeR+userRegionsSizeS);
	struct Quad quads[total_size];

	k=0;
	for(j=0;j<userRegionsSizeR;j++)   // Use x-dimension
	{
		struct Quad startTemp = {userRegionsR[j].xStart,j,user1,1};  // Add start of interval
		quads[k] = startTemp;
		k++;
		struct Quad endTemp = {userRegionsR[j].xEnd,j,user1,0};  // Add end of interval
		quads[k] = endTemp;
		k++;
	}
	for(j=0;j<userRegionsSizeS;j++)   // Use x-dimension
	{
		struct Quad startTemp = {userRegionsS[j].xStart,j,user2,1};  // Add start of interval
		quads[k] = startTemp;
		k++;
		struct Quad endTemp = {userRegionsS[j].xEnd,j,user2,0};  // Add end of interval
		quads[k] = endTemp;
		k++;
	}
	qsort(quads,k,sizeof(Quad),compareQuads);

	struct node *Dr = NULL;
	Dr = insert(Dr, 0, 0);


	struct node *Ds = NULL;
	Ds = insert(Ds, 0, 0);

	prev = quads[0].value;
	simn = 0;

	for(j=0;j<k;j++)
	{
		struct node *er = findMinimum(Dr);
		struct node *es = findMinimum(Ds);
		struct node *nextEs = findSuccessor(Ds, es);
		es = nextEs;
		

		while(er != NULL && es != NULL)
		{
			if(er->value < es->value)
			{
				struct node *erNext;
				struct node *esPrev, *esNext;
				esPrev = findPredecesor(Ds, es);
				erNext = findSuccessor(Dr, er);

				if(erNext != NULL)
					up = findMin(erNext,es);
				else
					up = es;
				if(esPrev != NULL)
				{
					simn += (up->value - er->value)*(quads[j].value - prev) * er->frequency*esPrev->frequency;
				}
				er = erNext;
			}
			else
			{
				struct node *erPrev = NULL;
				struct node *esNext = NULL;
				erPrev = findPredecesor(Dr, er);
				esNext = findSuccessor(Ds, es);

				if(esNext != NULL)
				{
					up = findMin(esNext,er);
				}
				else
					up = er;
				if(erPrev != NULL)
				{
					simn += (up->value - es->value)*(quads[j].value - prev) * es->frequency* erPrev->frequency;
				}
				es = esNext;
			}
		}

		if(quads[j].source == user1)
			Dr = updateDStructureQuads(Dr, quads[j], userRegionsR);
		else
			Ds = updateDStructureQuads(Ds, quads[j], userRegionsS);

		prev = quads[j].value;
	}
	return simn/(normrsq*normssq);
}


float computeSimilarityAndNorm(int userRegionsSizeR, Region *userRegionsR, int userRegionsSizeS, Region *userRegionsS, int user1, int user2)
{
	int i,j,k;
	float prev, prevs, prevr, simn;
	struct node *up;
	float normr=0.0, norms=0.0;

	int total_size = 2*(userRegionsSizeR+userRegionsSizeS);
	struct Quad quads[total_size];
	k=0;
	for(j=0;j<userRegionsSizeR;j++)   // Use x-dimension
	{
		struct Quad startTemp = {userRegionsR[j].xStart,j,user1,1};  // Add start of interval
		quads[k] = startTemp;
		k++;
		struct Quad endTemp = {userRegionsR[j].xEnd,j,user1,0};  // Add end of interval
		quads[k] = endTemp;
		k++;
	}
	for(j=0;j<userRegionsSizeS;j++)   // Use x-dimension
	{
		struct Quad startTemp = {userRegionsS[j].xStart,j,user2,1};  // Add start of interval
		quads[k] = startTemp;
		k++;
		struct Quad endTemp = {userRegionsS[j].xEnd,j,user2,0};  // Add end of interval
		quads[k] = endTemp;
		k++;
	}

	qsort(quads,k,sizeof(Quad),compareQuads);

	struct node *Dr = NULL;
	Dr = insert(Dr, 0, 0);

	struct node *Ds = NULL;
	Ds = insert(Ds, 0, 0);
	prev = quads[0].value;
	for(j=0;j<k;j++)
	{
		if(quads[j].source == user1)
		{
			prevr = quads[j].value;
			break;
		}
	}
	for(j=0;j<k;j++)
	{
		if(quads[j].source == user2)
		{
			prevs = quads[j].value;
			break;
		}
	}
	simn = 0;

	for(j=0;j<k;j++)
	{
		struct node *er = findMinimum(Dr);
		struct node *es = findMinimum(Ds);

		struct node *nextEs = findSuccessor(Ds, es);
		es = nextEs;

		while(er != NULL && es != NULL)
		{
			if(er->value < es->value)
			{
				struct node *erNext;

				struct node *esPrev, *esNext;
				esPrev = findPredecesor(Ds, es);
				erNext = findSuccessor(Dr, er);


				if(erNext != NULL)
					up = findMin(erNext,es);
				else
					up = es;
				if(esPrev != NULL)
				{
					if(esPrev->count < 0 )
						esPrev->count = 1;
					if(er->count < 0 )
						er->count = 1;
					simn += (up->value - er->value)*(quads[j].value - prev) * er->frequency* esPrev->frequency;
				}
				er = erNext;

			}
			else
			{
				struct node *erPrev = NULL;
				struct node *esNext = NULL;

				erPrev = findPredecesor(Dr, er);
				esNext = findSuccessor(Ds, es);
				
				if(esNext != NULL)
				{
					up = findMin(esNext,er);
				}
				else
					up = er;
				if(erPrev != NULL)
				{
					simn += (up->value - es->value)*(quads[j].value - prev) * es->frequency* erPrev->frequency;
				}
				es = esNext;
			}
		}

		if(quads[j].source == user1)
		{
			computeNsqInorder(Dr, &normr, prevr, quads[j].value,Dr);
			Dr = updateDStructureQuads(Dr, quads[j], userRegionsR);
			prevr = quads[j].value;
		}
		else
		{
			computeNsqInorder(Ds, &norms, prevs, quads[j].value, Ds);
			Ds = updateDStructureQuads(Ds, quads[j], userRegionsS);
			prevs = quads[j].value;
		}
		prev = quads[j].value;
	}
	return simn/(sqrt(normr)*sqrt(norms));
}

void computePairsInorder(struct node* node, int **pair, int *countPair)
{
    if (node == NULL)
        return;
 
    /* first recur on left child */
		computePairsInorder(node->left, pair, countPair);
		(*countPair)++;
		*pair[(*countPair)] = node->count;
		computePairsInorder(node->right, pair, countPair);
}

void internalLoop(Region rec, int recItem, Region *pivot, int pivotSize, float *area)
{
	int counter = 0;
	while(counter < pivotSize && rec.xEnd >= pivot[counter].xStart)
	{
		if(rec.yStart > pivot[counter].yEnd || rec.yEnd < pivot[counter].yStart || rec.xStart > pivot[counter].xEnd || rec.xEnd < pivot[counter].xStart)
		{
			counter++;
			continue;
		}
		float minX = max(rec.xStart,pivot[counter].xStart);
		float minY = max(rec.yStart,pivot[counter].yStart);
		float maxX = min(rec.xEnd,pivot[counter].xEnd);
		float maxY = min(rec.yEnd,pivot[counter].yEnd);
		(*area) += (maxX - minX)*(maxY - minY);	
		counter++;	

	}
	
}

float spatialSimilarity(int userRegionsSizeR, Region *userRegionsR, int userRegionsSizeS, Region *userRegionsS, float normrsq, float normssq)
{
	int i=0,j=0;
	int currentR = 0, currentS = 0;
	float area=0;


	while(currentR < userRegionsSizeR && currentS < userRegionsSizeS)
	{
		internalLoop(userRegionsS[currentS], currentS, userRegionsR, userRegionsSizeR, &area);
		currentS++;
	}
	return area/(normrsq*normssq);

}


int main(int argc, char **argv)
{

	FILE *f;
	int *userRegionsSize, *queryRegionsSize;
	Region **userRegions, **queryRegions;
	float *indexedUsersNsq, *queryUsersNsq;
	int indexedUsers, queryUsers, k_value;
	int i, j;
	clock_t t;
	double time_taken;
	int K = 0;
	clock_t startTime;


	char c;
  	int runProcessingMethod = -1;
	while ((c = getopt(argc, argv, "nsk:ph?")) != -1)
    {
        switch (c)
        {
            case 's':
                runProcessingMethod = SIMILARITY;
                break;
            case 'n':
                runProcessingMethod = SIMILARITY_AND_NORM;
                break;
            case 'p':
                runProcessingMethod = SPATIAL_SIMILARITY;
                break;
			case 'k':
                K = atoi(optarg);
                break;
            case '?':
            case 'h':
            default:
                // usage();
                
                break;
        }
    }


	cout << "User Region dataset : " << argv[optind] << endl;
  	cout << "Query Region dataset : " << argv[optind+1] << endl;

	f = fopen(argv[optind],"r");
	readData(f,&userRegionsSize, &userRegions, &indexedUsers);
	fclose(f);

    indexedUsersNsq = new float[indexedUsers];
	
	f = fopen(argv[optind+1], "r");

	readData(f, &queryRegionsSize, &queryRegions, &queryUsers);

	fclose(f);

    queryUsersNsq = new float[queryUsers];
    

    /* Compute Norms for Indexed Users */
	startTime = clock();
    for (int i = 0; i < indexedUsers; i++)
    {
        indexedUsersNsq[i] = computeNorm(userRegionsSize[i], userRegions[i]);
    }
	clock_t userNormTime = clock() - startTime;

    /* Compute Norms for Query Users */
	startTime = clock();
    for (int i = 0; i < queryUsers; i++)
    {
        queryUsersNsq[i] = computeNorm(queryRegionsSize[i], queryRegions[i]);
    }
	clock_t queryNormTime = clock() - startTime;

	switch (runProcessingMethod)
    {
        case SIMILARITY:
		{
			printf("Similarity - Algorithm 2\n");
			cerr << "UserNormTime\tqueryNormTime\tqueryTime\ttopKTime\n";
			
			for(int i = 0; i < queryUsers; i++)
			{
				startTime = clock();
				vector<pair<float, int>> topKUsers;
				for(j = 0; j < indexedUsers; j++) {
					float sim = computeSimilarity(userRegionsSize[j], userRegions[j], queryRegionsSize[i], queryRegions[i], indexedUsersNsq[j], queryUsersNsq[i],0,1);
					topKUsers.push_back(make_pair(sim, j));
				}
				clock_t queryTime = clock() - startTime;

				startTime = clock();
				sort(topKUsers.begin(), topKUsers.end(),greater<pair<float, int> >());
				clock_t topKTime = clock() - startTime;
			}
			break;
		}
		case SIMILARITY_AND_NORM:
		{
			printf("Similarity and norm - Algorithm 3\n");
			t = clock();
			vector<pair<float, int>> topKUsers;
			for(int i = 0; i < queryUsers; i++)
			{
				startTime = clock();
				for(j = 0; j < indexedUsers; j++) {

					float sim = computeSimilarityAndNorm(userRegionsSize[j], userRegions[j], queryRegionsSize[i], queryRegions[i],0,1);
					topKUsers.push_back(make_pair(sim, j));
				}
				clock_t queryTime = clock() - startTime;
				startTime = clock();
				// cout << "For the query User " << i << ":" << endl;
				sort(topKUsers.begin(), topKUsers.end(),greater<pair<float, int> >());
				clock_t topKTime = clock() - startTime;
			}
			break;
		}
		case SPATIAL_SIMILARITY:
		{
			t = clock(); 
			vector<pair<float, int>> topKUsers;
			for(int i = 0; i < queryUsers; i++) 
			{
				startTime = clock();
				for(j = 0; j < indexedUsers; j++) {
					float sim = spatialSimilarity(userRegionsSize[j], userRegions[j], queryRegionsSize[i], queryRegions[i], indexedUsersNsq[j], queryUsersNsq[i]);
					topKUsers.push_back(make_pair(sim, j));
				}
				clock_t queryTime = clock() - startTime;

				startTime = clock();
				sort(topKUsers.begin(), topKUsers.end(),greater<pair<float, int> >());
				clock_t topKTime = clock() - startTime;
			}
			break;
		}
	}

} 