#include <iostream>
#include <cstdio>
#include <cstdlib>
#include "RTree.h"
#include <cstring>
#include <unistd.h>
#include "bst.h"


#define ITERATIVE                                                       0
#define BATCH_SEARCH                                                    1
#define USER_CENTRIC                                                    2
#define USER_CENTRIC_SORTED                                             3
using namespace std;

#define DIMENSIONS 4

typedef int ValueType;

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

bool compareRegion(Region i1, Region i2)
{
    return (i1.xStart < i2.xStart);
}


int compare (const void * a, const void * b)
{

  Triple *tripleA = (Triple *)a;
  Triple *tripleB = (Triple *)b;

  if(tripleA->value > tripleB->value)
  	return 1;
  else if (tripleA->value < tripleB->value)
  	return -1;
  else
  	if(tripleA->type >= tripleB->type)
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


int readData(FILE *f, int **userRegionsSize, Region ***userRegions, int *numUsers)
{
	char *line = NULL;
	size_t len = 0; 
	ssize_t read; 
	const char delim[2] = " ";
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

void internalLoop_sweepX(Region rec, Region *pivot, int pivotStart,int pivotSize, float *area)
{
	int counter = pivotStart;
	while(counter < pivotSize && rec.xEnd >= pivot[counter].xStart)
	{
		if(rec.yStart > pivot[counter].yEnd || rec.yEnd < pivot[counter].yStart)
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

void internalLoop(Region rec, Region *pivot, int pivotSize, float *area)
{
	int counter = 0;
	while(counter < pivotSize)
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

float spatialSimilarity(int userRegionsSizeR, Region *userRegionsR, int userRegionsSizeS, Region *userRegionsS, float normrsq, float normssq, int flag)
{
	int i=0,j=0;
	int currentR = 0, currentS = 0;
	float area=0;

    switch (flag)
    {
    case USER_CENTRIC:
        while(currentR < userRegionsSizeR)
	    {
            internalLoop(userRegionsR[currentR],userRegionsS, userRegionsSizeS, &area);
            currentR++;
        }
    
        break;
    case USER_CENTRIC_SORTED:    
        while(currentR < userRegionsSizeR && currentS < userRegionsSizeS)
        {
            if(currentR < currentS)
            {
                internalLoop_sweepX(userRegionsR[currentR], userRegionsS, currentS,userRegionsSizeS, &area);
                currentR++;
            }
            else
            {
                internalLoop_sweepX(userRegionsS[currentS], userRegionsR, currentR,userRegionsSizeR, &area);
                currentS++;
            }
            
        }
        break;
    }
	return area/(normrsq*normssq);

}

int main(int argc, char **argv)
{
    FILE *f;
	int *userRegionsSize, *queryRegionsSize;
	Region **userRegions, **queryRegions;
	float *indexedUsersNsq, *queryUsersNsq;
    int maxNumRegions=-1;
	int indexedUsers, queryUsers, k_value;
    char c;
    int indexMethod = -1, numIterations = 0;
	clock_t startTime, indexingTime;
	double time_taken;

    typedef RTree<ValueType, float, 2, float> MyTree;
    
    MyTree tree;
    
    while ((c = getopt(argc, argv, "n:b:u:s:k:i:")) != -1)
    {
        switch(c)
        {
          
            case 'n':
                indexMethod = ITERATIVE;
                break;
            case 'b':
                indexMethod = BATCH_SEARCH;
                break;
            case 'u':
                indexMethod = USER_CENTRIC;
                break;
            case 's':
                indexMethod = USER_CENTRIC_SORTED;
                break;
            case 'k':
                k_value = atoi(optarg);
                break;
            case 'i':
                numIterations = atoi(optarg);
                break;
            default:
                cout << "Invalid arguments!" << endl;
            break;
        }
    }

    cout << "Indexing Method\t" << "Indexed Users\t" << "Query Users"<< endl;
    cout << indexMethod << "\t\t" << argv[optind-1] << "\t" << argv[optind] << endl;
    
    f = fopen(argv[optind-1],"r");
	readData(f,&userRegionsSize, &userRegions, &indexedUsers);
	fclose(f);

	printf("Total Indexed users = %d\n",indexedUsers);
    indexedUsersNsq = new float[indexedUsers];
	
	f = fopen(argv[optind], "r");

	readData(f, &queryRegionsSize, &queryRegions, &queryUsers);

	fclose(f);
	printf("Total Query users = %d\n\n", queryUsers);
    
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
    
    // cout << "NSQ = " << indexedUsersNsq[0] << endl;

    switch (indexMethod)
    {
        case ITERATIVE:
            
            cerr << "UserNormTime\tqueryNormTime\tindexingTime\tqueryTime\ttopKTime\n";
            startTime = clock();
            initialize_user_metrics(indexedUsersNsq, indexedUsers);

            for(int i = 0; i < indexedUsers; i++) {
                for (int j = 0; j < userRegionsSize[i]; j++)
                {
                    Rect r = Rect(userRegions[i][j].xStart, userRegions[i][j].yStart, userRegions[i][j].xEnd, userRegions[i][j].yEnd);
                    tree.Insert(r.min, r.max, i);
                }
            }
            indexingTime = clock() - startTime;

            
            for (int i = 0; i < queryUsers; i++)
            {
                startTime = clock();


                for (int j = 0; j < queryRegionsSize[i]; j++)
                {   
                    Rect q = Rect(queryRegions[i][j].xStart, queryRegions[i][j].yStart, queryRegions[i][j].xEnd, queryRegions[i][j].yEnd);
                    tree.Search(q.min, q.max);
                }
                clock_t queryTime = clock() - startTime;


                startTime = clock();
                vector<pair<float, int>> topKUsers = top_K_Users(indexedUsers, k_value, queryUsersNsq[i]);
                clock_t topKTime = clock() - startTime;
                
                cerr << ((double)userNormTime)/CLOCKS_PER_SEC << "\t\t" << ((double)queryNormTime)/CLOCKS_PER_SEC << "\t\t" << ((double)indexingTime)/CLOCKS_PER_SEC << "\t\t" << ((double)queryTime)/CLOCKS_PER_SEC << "\t\t" << ((double)topKTime)/CLOCKS_PER_SEC << endl;

                cout << "For the query User " << i << ":" << endl;
                for (auto it= topKUsers.begin(); it != topKUsers.end(); it++)
                {
                    cout << "\tUser ID: " << (*it).second << " - similarity: " << (*it).first << endl; 
                }
                cout << "---------------------------------------------------------------------" << endl<< endl;    
            }

            break;
        
        case BATCH_SEARCH:

            cerr << "UserNormTime\tqueryNormTime\tindexingTime\tqueryTime\ttopKTime\n";

            startTime = clock();
            initialize_user_metrics(indexedUsersNsq, indexedUsers);

            for(int i = 0; i < indexedUsers; i++) {
                for (int j = 0; j < userRegionsSize[i]; j++)
                {
                    Rect r = Rect(userRegions[i][j].xStart, userRegions[i][j].yStart, userRegions[i][j].xEnd, userRegions[i][j].yEnd);
                    tree.Insert(r.min, r.max, i);
                }
            }
            indexingTime = clock() - startTime;

            float minXQueryMBR, minYQueryMBR;
            float maxXQueryMBR, maxYQueryMBR;

            for (int i =0; i < queryUsers; i++)
            {
                startTime = clock();

                minXQueryMBR = queryRegions[i][0].xStart;
                minYQueryMBR = queryRegions[i][0].yStart;
                maxXQueryMBR = queryRegions[i][0].xEnd;
                maxYQueryMBR = queryRegions[i][0].yEnd;

                for (int j = 1; j < queryRegionsSize[i]; j++)
                {
                    minXQueryMBR = min(minXQueryMBR,queryRegions[i][j].xStart);
                    minYQueryMBR = min(minYQueryMBR, queryRegions[i][j].yStart);
                    maxXQueryMBR = max(maxXQueryMBR, queryRegions[i][j].xEnd);
                    maxYQueryMBR = max(maxYQueryMBR, queryRegions[i][j].yEnd);
                }

                Rect q = Rect(minXQueryMBR, minYQueryMBR, maxXQueryMBR, maxYQueryMBR);
                tree.Batch_Search(q.min, q.max, queryRegions[i], queryRegionsSize[i]);

                clock_t queryTime = clock() - startTime;

                startTime = clock();
                vector<pair<float, int>> topKUsers = top_K_Users(indexedUsers, k_value, queryUsersNsq[i]);
                clock_t topKTime = clock() - startTime;

                cerr << ((double)userNormTime)/CLOCKS_PER_SEC << "\t\t" << ((double)queryNormTime)/CLOCKS_PER_SEC << "\t\t" << ((double)indexingTime)/CLOCKS_PER_SEC << "\t\t" << ((double)queryTime)/CLOCKS_PER_SEC << "\t\t" << ((double)topKTime)/CLOCKS_PER_SEC << endl;

            
                cout << "For the query User " << i << ":" << endl;
                for (auto it= topKUsers.begin(); it != topKUsers.end(); it++)
                {
                    cout << "\tUser ID: " << (*it).second << " - similarity: " << (*it).first << endl; 
                }
                cout << "---------------------------------------------------------------------" << endl<< endl;
            } 
            break;
        case USER_CENTRIC:

            cerr << "UserNormTime\tqueryNormTime\tindexingTime\tqueryTime\ttopKTime\n";

            startTime = clock();
            for(int i = 0; i < indexedUsers; i++) {
        
                float minXUserMBR, minYUserMBR;
                float maxXUserMBR, maxYUserMBR;

                minXUserMBR = userRegions[i][0].xStart;
                minYUserMBR = userRegions[i][1].yStart;
                maxXUserMBR = userRegions[i][2].xEnd;
                maxYUserMBR = userRegions[i][3].yEnd;


                for (int j = 1; j < userRegionsSize[i]; j++)
                {
                    minXUserMBR = min(minXUserMBR, userRegions[i][j].xStart);
                    minYUserMBR = min(minYUserMBR, userRegions[i][j].yStart);
                    maxXUserMBR = max(maxXUserMBR, userRegions[i][j].xEnd);
                    maxYUserMBR = max(maxYUserMBR, userRegions[i][j].yEnd);
                }

                Rect r = Rect(minXUserMBR, minYUserMBR, maxXUserMBR, maxYUserMBR);
                tree.Insert(r.min, r.max,i);
            }

            indexingTime = clock() - startTime;

            for (int i =0; i < queryUsers; i++)
            {
                startTime = clock();

                minXQueryMBR = queryRegions[i][0].xStart;
                minYQueryMBR = queryRegions[i][1].yStart;
                maxXQueryMBR = queryRegions[i][2].xEnd;
                maxYQueryMBR = queryRegions[i][3].yEnd;

                for (int j = 1; j < queryRegionsSize[i]; j++)
                {
                    minXQueryMBR = min(minXQueryMBR,queryRegions[i][j].xStart);
                    minYQueryMBR = min(minYQueryMBR, queryRegions[i][j].yStart);
                    maxXQueryMBR = max(maxXQueryMBR, queryRegions[i][j].xEnd);
                    maxYQueryMBR = max(maxYQueryMBR, queryRegions[i][j].yEnd);
                }

                Rect q = Rect(minXQueryMBR, minYQueryMBR, maxXQueryMBR, maxYQueryMBR);

                vector<int> hits = tree.Batch_Search2(q.min, q.max, queryRegions[i], queryRegionsSize[i]);


                cout << "For the query User " << i << ":" << endl;
                if (hits.size() <= k_value)
                {
                    for(int j = 0; j < hits.size(); j++)
                    {
                        int idx = hits[j];
                        float similarity = spatialSimilarity(queryRegionsSize[i], queryRegions[i], userRegionsSize[idx], userRegions[idx], queryUsersNsq[i], indexedUsersNsq[idx], USER_CENTRIC);
                        clock_t queryTime = clock() - startTime;
                        cerr << ((double)userNormTime)/CLOCKS_PER_SEC << "\t\t" << ((double)queryNormTime)/CLOCKS_PER_SEC << "\t\t" << ((double)indexingTime)/CLOCKS_PER_SEC << "\t\t" << ((double)queryTime)/CLOCKS_PER_SEC << "\t\t0" << endl;
                    }
                }
                else
                {
                    vector<pair<float, int>> topKUsers;
                    for(int j = 0; j < hits.size(); j++)
                    {
                        int idx = hits[j];
                        float similarity = spatialSimilarity(queryRegionsSize[i], queryRegions[i], userRegionsSize[idx], userRegions[idx], queryUsersNsq[i], indexedUsersNsq[idx], USER_CENTRIC);
                        topKUsers.push_back(make_pair(similarity, idx));
                    }

                    sort(topKUsers.begin(), topKUsers.end(),sortbysec);
                    clock_t queryTime = clock() - startTime;
                    cerr << ((double)userNormTime)/CLOCKS_PER_SEC << "\t\t" << ((double)queryNormTime)/CLOCKS_PER_SEC << "\t\t" << ((double)indexingTime)/CLOCKS_PER_SEC << "\t\t" << ((double)queryTime)/CLOCKS_PER_SEC << "\t\t0" << endl;
 
                    for (int j = 0; j < k_value; j++)
                    {
                        cout << "\tUser ID: " << topKUsers[j].second << " - similarity: " << topKUsers[j].first << endl;
                    }
                }

                cout << "---------------------------------------------------------------------" << endl<< endl;
            }
            break;
        case USER_CENTRIC_SORTED:
            cerr << "UserNormTime\tqueryNormTime\tindexingTime\tqueryTime\ttopKTime\n";


            startTime = clock();
            for(int i = 0; i < indexedUsers; i++) {
        
                float minXUserMBR, minYUserMBR;
                float maxXUserMBR, maxYUserMBR;

                sort(userRegions[i], userRegions[i]+userRegionsSize[i], compareRegion);
                
                minXUserMBR = userRegions[i][0].xStart;
                minYUserMBR = userRegions[i][1].yStart;
                maxXUserMBR = userRegions[i][2].xEnd;
                maxYUserMBR = userRegions[i][3].yEnd;

                for (int j = 1; j < userRegionsSize[i]; j++)
                {
                    minXUserMBR = min(minXUserMBR, userRegions[i][j].xStart);
                    minYUserMBR = min(minYUserMBR, userRegions[i][j].yStart);
                    maxXUserMBR = max(maxXUserMBR, userRegions[i][j].xEnd);
                    maxYUserMBR = max(maxYUserMBR, userRegions[i][j].yEnd);
                }

                Rect r = Rect(minXUserMBR, minYUserMBR, maxXUserMBR, maxYUserMBR);
                tree.Insert(r.min, r.max,i);
            }
            indexingTime = clock() - startTime;

            for (int i =0; i < queryUsers; i++)
            {
                startTime = clock();
                sort(queryRegions[i], queryRegions[i]+queryRegionsSize[i], compareRegion);


                minXQueryMBR = queryRegions[i][0].xStart;
                minYQueryMBR = queryRegions[i][1].yStart;
                maxXQueryMBR = queryRegions[i][2].xEnd;
                maxYQueryMBR = queryRegions[i][3].yEnd;

                for (int j = 1; j < queryRegionsSize[i]; j++)
                {
                    minXQueryMBR = min(minXQueryMBR,queryRegions[i][j].xStart);
                    minYQueryMBR = min(minYQueryMBR, queryRegions[i][j].yStart);
                    maxXQueryMBR = max(maxXQueryMBR, queryRegions[i][j].xEnd);
                    maxYQueryMBR = max(maxYQueryMBR, queryRegions[i][j].yEnd);
                }

                Rect q = Rect(minXQueryMBR, minYQueryMBR, maxXQueryMBR, maxYQueryMBR);

                vector<int> hits = tree.Batch_Search2(q.min, q.max, queryRegions[i], queryRegionsSize[i]);

                cout << "For the query User " << i << ":" << endl;
                if (hits.size() <= k_value)
                {   
                    for(int j = 0; j < hits.size(); j++)
                    {
                        int idx = hits[j];
                        float similarity = spatialSimilarity(queryRegionsSize[i], queryRegions[i], userRegionsSize[idx], userRegions[idx], queryUsersNsq[i], indexedUsersNsq[idx], USER_CENTRIC_SORTED);
                        clock_t queryTime = clock() - startTime;
                        cerr << ((double)userNormTime)/CLOCKS_PER_SEC << "\t\t" << ((double)queryNormTime)/CLOCKS_PER_SEC << "\t\t" << ((double)indexingTime)/CLOCKS_PER_SEC << "\t\t" << ((double)queryTime)/CLOCKS_PER_SEC << "\t\t0" << endl;
                    }
                }
                else
                {
                    vector<pair<float, int>> topKUsers;
                    for(int j = 0; j < hits.size(); j++)
                    {
                        int idx = hits[j];
                        float similarity = spatialSimilarity(queryRegionsSize[i],queryRegions[i], userRegionsSize[idx], userRegions[idx],queryUsersNsq[i], indexedUsersNsq[idx], USER_CENTRIC_SORTED);
                        topKUsers.push_back(make_pair(similarity, idx));
                    }

                    sort(topKUsers.begin(), topKUsers.end(),sortbysec);
                    
                    clock_t queryTime = clock() - startTime;
                    cerr << ((double)userNormTime)/CLOCKS_PER_SEC << "\t\t" << ((double)queryNormTime)/CLOCKS_PER_SEC << "\t\t" << ((double)indexingTime)/CLOCKS_PER_SEC << "\t\t" << ((double)queryTime)/CLOCKS_PER_SEC << "\t\t0" << endl;
                    
                    for (int j = 0; j < k_value; j++)
                    {
                        cout << "\tUser ID: " << topKUsers[j].second << " - similarity: " << topKUsers[j].first << endl;
                    }
                }
               

                cout << "---------------------------------------------------------------------" << endl<< endl;
            }
            break;
    }

}