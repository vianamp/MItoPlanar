#include <list>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include <igraph/igraph.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkCellArray.h>
#include <vtkBitArray.h>
#include <vtkKdTreePointLocator.h>

#define _min(a,b) ((a<b) ? a : b)
#define _max(a,b) ((a>b) ? a : b)
#define _eps 1E-6
#define TWO_PI 6.2831853071795864769252866

/* =================================================================
   CLASSES
   =================================================================*/

struct _clock{
	long long int t, nreset;
	_clock() {
		t = nreset = 0;
	}
	void reset() {
		t = nreset = 0;	
	}
};

struct _node{
	double x;
	double y;
};

struct _edge{
	int i;
	int j;
	int type;
	double length;
	bool operator == (const _edge &edge) {
		return ((i==edge.i)&&(j==edge.j))||((i==edge.j)&&(j==edge.i));
	};
};

class _Graph{
private:
public:
	int N;
	float Lx, Ly;
	std::list<_edge> Edges;
	std::vector<_node> Nodes;
	_clock clock;

	_Graph() {
		N = 0;
		Lx = Ly = 0.0;
	};

	void SavePolyData(const char FileName[]);
	void SaveGNET(const char FileName[]);
	void MakeShallowCopy(const char FilePrefix[], int *NEdges, double *Length, double *pk1);
	void MakeDeepCopy(const char FilePrefix[], int *NEdges, double *Length, double *pk1);
	void GenerateFromScratch(int n, int sn);

	double GetEdgeLength(int i, int j);
	bool IsPlanarEdge(int i, int j);
	int GetTotalNumberOfCrosses();
	void CreateVirtualNodes();
	int GetClosestVirtualNode(int i, int j);
	void ClipNetwork();
	void DeleteEdge(_edge edge);
	void DeleteRandomEdge(int *, int *j);
	void GetProperties(double *M);
	void GetMinimumDistanceVector(std::vector<double> &Dmin);
	void ShuffleCoordinates();
	bool SwapPairOfEdges(int *i, int *j, int *p, int *q);
};

/* =================================================================
   AUXILIAR FUNCTIONS
   =================================================================*/

void _help() {
	printf("-path\n");
	printf("-model\t\t{EDG,LGT,PK1,WAX,WXS}\n");
	printf("-r\t\t(number of instances)\n");
	printf("-alpha\t\t(distance strength in distance-based models)\n");
	printf("-save\t\t(save instances generated)\n");
	printf("-pk1_mode\n");
	printf("-size_mode\n");
}

double generateGaussianNoise(const double &variance) {
	static bool haveSpare = false;
	static double rand1, rand2;
 	if(haveSpare) {
		haveSpare = false;
		return sqrt(variance * rand1) * sin(rand2);
	}
 
	haveSpare = true;
 	rand1 = rand() / ((double) RAND_MAX);
	if(rand1 < 1e-100) rand1 = 1e-100;
	rand1 = -2 * log(rand1);
	rand2 = (rand() / ((double) RAND_MAX)) * TWO_PI;
 	return sqrt(variance * rand1) * cos(rand2);
}

int _mirror(int q) {
	switch (q) {
		case 1: return 3;
		case 2: return 4;
		case 3: return 1;
	 default  : return 2;
	}
}

bool IsInBetween(long double ri[3], long double rj[3], long double v, int d) {
	return ( v > _min(ri[d],rj[d])+_eps ) && ( v < _max(ri[d],rj[d]-_eps) ) ? true : false;
}

bool BoundBoxIntercept(long double ri[3], long double rj[3], long double ru[3], long double rv[3]) {
	if (_min(ri[0],rj[0])<_max(ru[0],rv[0])) {
		if (_max(ri[0],rj[0])>_min(ru[0],rv[0])) {
			if (_min(ri[1],rj[1])<_max(ru[1],rv[1])) {
				if (_max(ri[1],rj[1])>_min(ru[1],rv[1])) return true;
			}
		}
	}
	return false;
}

bool SegmentsContainPoint(long double ri[3], long double rj[3], long double ru[3], long double rv[3], long double xc, long double yc) {
	if ( IsInBetween(ri,rj,xc,0) ) {
		if ( IsInBetween(ri,rj,yc,1) ) {
			if ( IsInBetween(ru,rv,xc,0) ) {
				if ( IsInBetween(ru,rv,yc,1) ) {
					return true;
				}
			}
		}
	}
	return false;
}

int _rng(int i) {
	return rand() % i;
}

/* =================================================================
   I/O METHODS
   =================================================================*/

void _Graph::SavePolyData(const char FileName[]) {
	printf("Saving polydata...\n");

	vtkSmartPointer<vtkPoints> Points = vtkSmartPointer<vtkPoints>::New();
	for (int i=0; i < Nodes.size(); i++) {
		Points -> InsertNextPoint(Nodes[i].x,Nodes[i].y,0.0);
	}

	vtkSmartPointer<vtkCellArray> Cells = vtkSmartPointer<vtkCellArray>::New();
	for (std::list<_edge>::iterator it = Edges.begin(), end = Edges.end(); it != end; ++it) {
		Cells -> InsertNextCell(2);
		Cells -> InsertCellPoint((*it).i);
		Cells -> InsertCellPoint((*it).j);
	}

	vtkSmartPointer<vtkPolyData> Net = vtkSmartPointer<vtkPolyData>::New();
	Net -> SetPoints(Points);
	Net -> SetLines(Cells);

	vtkSmartPointer<vtkPolyDataWriter> PolyWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
	PolyWriter -> SetFileName(FileName);
	PolyWriter -> SetInputData(Net);
	PolyWriter -> Write();
	printf("\tDone!\n");
}

void _Graph::SaveGNET(const char FileName[]) {
	_edge edge;
	FILE *f = fopen(FileName,"w");
	fprintf(f,"%d\n",N);
	for (std::list<_edge>::iterator it = Edges.begin(), end = Edges.end(); it != end; ++it) {
		edge = *it;
		fprintf(f,"%d\t%d\t%1.6f\n",edge.i,edge.j,edge.length);
	}
	fclose(f);
}

void _Graph::MakeShallowCopy(const char FilePrefix[], int *NEdges, double *Length, double *pk1) {

	#ifdef DEBUG
		printf("Copying Coordinates From %s.coo2d\n",FilePrefix);
	#endif

	char text[128], _path[128];
	int i, j, E = 0;
	float x, y, length, L = 0.0;
	sprintf(_path,"%s.gnet",FilePrefix);
	FILE *fg = fopen(_path,"r");
	if (!fg) {
		printf("File %s.gnet not found.",FilePrefix);
	}
	sprintf(_path,"%s.coo2d",FilePrefix);
	FILE *fc = fopen(_path,"r");
	if (!fc) {
		printf("File %s.coo2d not found.",FilePrefix);
	}
	fscanf(fg,"%d",&N);
	std::vector<int> K (N,0);
	while (fscanf(fg,"%d %d %f",&i,&j,&length)!=EOF) {
		K[i]++;
		K[j]++;
		E++;
		L+=length;
	}
	fclose(fg);

	Nodes.clear();

	#ifdef DEBUG
		printf("\t#Nodes = %d\n",N);
		printf("\t#Edges = %d\n",E);
		printf("\t#Length = %1.3f\n",L);
	#endif

	fscanf(fc,"%s",text);
	fscanf(fc,"%f",&Lx);
	fscanf(fc,"%f",&Ly);
	fscanf(fc,"%s",text);

	while (fscanf(fc,"%f %f",&x,&y)!=EOF) {
		x += _eps * (1.0*rand())/RAND_MAX;
		y += _eps * (1.0*rand())/RAND_MAX;
		_node node = {x,y};
		Nodes.push_back(node);
	}
	fclose(fc);
	*Length = L;
	*NEdges = E;

	int nk1 = 0;
	for (i = 0; i < N; i++) if (K[i]==1) nk1++;

	*pk1 = (double)nk1 / N;
	K.clear();

	#ifdef DEBUG
		printf("\tLx = %1.3f, Ly = %1.3f\n",Lx,Ly);
		printf("\tCopy Complete!\n");
	#endif
}

void _Graph::MakeDeepCopy(const char FilePrefix[], int *NEdges, double *Length, double *pk1) {

	#ifdef DEBUG
		printf("Copying Coordinates From %s.coo2d\n",FilePrefix);
	#endif

	char text[128], _path[128];
	int i, j, E = 0;
	float x, y, length, L = 0.0;
	sprintf(_path,"%s.gnet",FilePrefix);
	FILE *fg = fopen(_path,"r");
	if (!fg) {
		printf("File %s.gnet not found.",FilePrefix);
	}
	sprintf(_path,"%s.coo2d",FilePrefix);
	FILE *fc = fopen(_path,"r");
	if (!fc) {
		printf("File %s.coo2d not found.",FilePrefix);
	}
	fscanf(fg,"%d",&N);
	std::vector<int> K (N,0);
	while (fscanf(fg,"%d %d %f",&i,&j,&length)!=EOF) {
		K[i]++;
		K[j]++;
		E++;
		L+=length;
	}

	Nodes.clear();

	#ifdef DEBUG
		printf("\t#Nodes = %d\n",N);
		printf("\t#Edges = %d\n",E);
		printf("\t#Length = %1.3f\n",L);
	#endif

	fscanf(fc,"%s",text);
	fscanf(fc,"%f",&Lx);
	fscanf(fc,"%f",&Ly);
	fscanf(fc,"%s",text);

	while (fscanf(fc,"%f %f",&x,&y)!=EOF) {
		x += _eps * (1.0*rand())/RAND_MAX;
		y += _eps * (1.0*rand())/RAND_MAX;
		_node node = {x,y};
		Nodes.push_back(node);
	}
	fclose(fc);
	*Length = L;
	*NEdges = E;

	int nk1 = 0;
	for (i = 0; i < N; i++) if (K[i]==1) nk1++;

	*pk1 = (double)nk1 / N;
	K.clear();

	#ifdef DEBUG
		printf("\tLx = %1.3f, Ly = %1.3f\n",Lx,Ly);
	#endif

	CreateVirtualNodes();

	long int k;
	int q, ne = 0;

	#ifdef DEBUG
		printf("Copying edges form %s\n",FilePrefix);
	#endif

	Edges.clear();

	rewind(fg);
	fscanf(fg,"%d",&i);

	int _nplanar=0;
	while (fscanf(fg,"%d %d %f",&i,&j,&length) != EOF) {
		q = GetClosestVirtualNode(i,j);
		length = GetEdgeLength(i,j + q*N);
		if (q) {
			_edge edge = {i,j+q*N,5,length};
			Edges.push_back(edge);
			_edge edge_m = {i+_mirror(q)*N,j,6,length};
			Edges.push_back(edge_m);
		} else {
			for (q = 0; q < 5; q++) {
	 			_edge edge = {i+q*N,j+q*N,q,length};
	 			Edges.push_back(edge);
	 		}
		}
	}

	fclose(fg);

}

void _Graph::GenerateFromScratch(int n, int sn) {
	#ifdef DEBUG
		printf("Generating from scratch...\n");
	#endif

	N = 0;
	while (N < 13)
		N = n + (int)generateGaussianNoise((double)sn);

	Lx = 15;

	Ly = 8;

	Nodes.clear();

	#ifdef DEBUG
		printf("\t#Nodes = %d\n",N);
	#endif

	for (int i = 0; i < N; i++) {
		_node node = {0,0};
		Nodes.push_back(node);
	}
	
	ShuffleCoordinates();

	#ifdef DEBUG
		printf("\tLx = %1.3f, Ly = %1.3f\n",Lx,Ly);
		printf("\tComplete!\n");
	#endif

}

/* =================================================================
   AUXILIAR METHODS
   =================================================================*/

double _Graph::GetEdgeLength(int i, int j) {
	double ri[3], rj[3];
	ri[0] = Nodes[i].x;
	ri[1] = Nodes[i].y;
	rj[0] = Nodes[j].x;
	rj[1] = Nodes[j].y;
	return sqrt(pow(ri[0]-rj[0],2)+pow(ri[1]-rj[1],2));
}

bool _Graph::IsPlanarEdge(int i, int j) {
	long double Det, Q, xc, yc, fac;
	long double ri[3], rj[3], ru[3], rv[3];
	ri[0] = Nodes[i].x; ri[1] = Nodes[i].y;
	rj[0] = Nodes[j].x; rj[1] = Nodes[j].y;

	bool _stillplanar = true;

	#pragma omp parallel
	for (std::list<_edge>::iterator it = Edges.begin(), end = Edges.end(); (it!=end)&&(_stillplanar); ++it) {
		ru[0] = Nodes[(*it).i].x; ru[1] = Nodes[(*it).i].y;
		rv[0] = Nodes[(*it).j].x; rv[1] = Nodes[(*it).j].y;

		if (BoundBoxIntercept(ri,rj,ru,rv)) {

			Det = ( rj[1] - ri[1] ) * ( rv[0] - ru[0] ) - ( rv[1] - ru[1] ) * ( rj[0] - ri[0] );

			if(fabs(Det) > _eps) {
				Q =  ri[0] * ( rv[1] - ru[1] ) + ru[0] * ( ri[1] - rv[1] ) + rv[0] * ( ru[1] - ri[1] );
				xc = ri[0] + (Q / Det) * ( rj[0] - ri[0] );
				yc = ri[1] + (Q / Det) * ( rj[1] - ri[1] );

				//if (SegmentsContainPoint(ri,rj,ru,rv,xc,yc)) return false;
				if (SegmentsContainPoint(ri,rj,ru,rv,xc,yc)) _stillplanar = false;
			}

		}
	}
	return _stillplanar;
}

int _Graph::GetTotalNumberOfCrosses() {
	int _ncrosses = 0;
	long double Det, Q, xc, yc, fac;
	long double ri[3], rj[3], ru[3], rv[3];

	#pragma omp parallel
	for (std::list<_edge>::iterator it_o = Edges.begin(), end_o = Edges.end(); (it_o!=end_o); ++it_o) {

		if (  ((*it_o).i<N && (*it_o).j<N) || (((*it_o).i<N || (*it_o).i<N) && ((*it_o).i<(*it_o).j))  ) {

			ri[0] = Nodes[(*it_o).i].x; ri[1] = Nodes[(*it_o).i].y;
			rj[0] = Nodes[(*it_o).j].x; rj[1] = Nodes[(*it_o).j].y;

			for (std::list<_edge>::iterator it = Edges.begin(), end = Edges.end(); (it!=end); ++it) {

				if ((*it_o).i!=(*it).i||(*it_o).j!=(*it).j) {

					ru[0] = Nodes[(*it).i].x; ru[1] = Nodes[(*it).i].y;
					rv[0] = Nodes[(*it).j].x; rv[1] = Nodes[(*it).j].y;

					if (BoundBoxIntercept(ri,rj,ru,rv)) {

						Det = ( rj[1] - ri[1] ) * ( rv[0] - ru[0] ) - ( rv[1] - ru[1] ) * ( rj[0] - ri[0] );

						if(fabs(Det) > _eps) {
							Q =  ri[0] * ( rv[1] - ru[1] ) + ru[0] * ( ri[1] - rv[1] ) + rv[0] * ( ru[1] - ri[1] );
							xc = ri[0] + (Q / Det) * ( rj[0] - ri[0] );
							yc = ri[1] + (Q / Det) * ( rj[1] - ri[1] );

							if (SegmentsContainPoint(ri,rj,ru,rv,xc,yc)) _ncrosses++;

						}

					}

				}

			}

		}
	}
	return _ncrosses/2;
}

void _Graph::CreateVirtualNodes() {
	double r[3], x, y;
	N = Nodes.size();
	double Q[4][2] = {{0,Ly},{Lx,0},{0,-Ly},{-Lx,0}};
	for (int q = 1; q <= 4; q++) {
		for (int i = 0; i < N; i++) {
			_node node = {Nodes[i].x+Q[q-1][0],Nodes[i].y+Q[q-1][1]};
			Nodes.push_back(node);
		}
	}
}

int _Graph::GetClosestVirtualNode(int i, int j) {
	int qo;
	double ri[3], rj[3], d, dmin = 1E5;
	ri[0] = Nodes[i].x;
	ri[1] = Nodes[i].y;
	for (int q = 0; q < 5; q++) {
		rj[0] = Nodes[j+q*N].x;
		rj[1] = Nodes[j+q*N].y;
		d = pow(rj[1]-ri[1],2) + pow(rj[0]-ri[0],2);
		if (d < dmin) {
			dmin = d;
			qo = q;
		}
	}
	return qo;
}

void _Graph::ClipNetwork() {

	#ifdef DEBUG
		printf("Clipping network...\n");
	#endif

	int i, j;
	double d;
	_edge edge;
	for (std::list<_edge>::iterator it = Edges.begin(), end = Edges.end(); it != end; ++it) {
		edge = *it;
		i = edge.i;
		j = edge.j;
		if ( edge.type ) {
			if ( edge.type == 5 ) {
				(*it).j = j-(j/N)*N;
				(*it).type = 0;
				(*it).length = sqrt(pow(Nodes[i].x-Nodes[j].x,2)+pow(Nodes[i].y-Nodes[j].y,2));
			} else Edges.erase(it--);
		} else {
			d = sqrt(pow(Nodes[i].x-Nodes[j].x,2)+pow(Nodes[i].y-Nodes[j].y,2));
			(*it).length = d;
		}
	}

	for (i=N;i<5*N;i++) Nodes.pop_back();

	#ifdef DEBUG
		printf("\tDone!\n");
	#endif
}

void _Graph::DeleteEdge(_edge edge) {

	int q, i, j;
	std::list<_edge>::iterator it;

	i = edge.i;
	j = edge.j;

	if (edge.type < 5) {

		i -= (i/N)*N;
		j -= (j/N)*N;

		for (q = 0; q < 5; q++) {
			edge.i = i + q*N;
			edge.j = j + q*N;
			it = std::find(Edges.begin(),Edges.end(),edge);
			Edges.erase(it);
		}

	} else if (edge.type == 5) {

		q = j / N;

		Edges.erase(it);
		edge.i = i+_mirror(j/N)*N;
		edge.j = j - (j/N)*N;
		it = std::find(Edges.begin(),Edges.end(),edge);
		Edges.erase(it);

	} else {

		Edges.erase(it);
		edge.i = i - (i/N)*N;
 		edge.j = j + _mirror(i/N)*N;
		it = std::find(Edges.begin(),Edges.end(),edge);
		Edges.erase(it);

	}

}

void _Graph::DeleteRandomEdge(int *io, int *jo) {

	int q, i, j;
	std::list<_edge>::iterator it = Edges.begin();
	std::advance(it,rand()%Edges.size());
	_edge edge = *it;

	i = edge.i;
	j = edge.j;

	if (edge.type < 5) {

		i -= (i/N)*N;
		j -= (j/N)*N;

		*io = i - (i/N)*N;
		*jo = j - (j/N)*N;

		for (q = 0; q < 5; q++) {
			edge.i = i + q*N;
			edge.j = j + q*N;
			it = std::find(Edges.begin(),Edges.end(),edge);
			Edges.erase(it);
		}

	} else if (edge.type == 5) {

		q = j / N;

		*io = i;
		*jo = j - (j/N)*N;
		Edges.erase(it);
		edge.i = i+_mirror(j/N)*N;
		edge.j = j - (j/N)*N;
		it = std::find(Edges.begin(),Edges.end(),edge);
		Edges.erase(it);

	} else {

		*io = i - (i/N)*N;
		*jo = j;
		Edges.erase(it);
		edge.i = i - (i/N)*N;
 		edge.j = j + _mirror(i/N)*N;
		it = std::find(Edges.begin(),Edges.end(),edge);
		Edges.erase(it);

	}

}

void _Graph::GetProperties(double *M) {

	igraph_i_set_attribute_table(&igraph_cattribute_table);

	igraph_t iGraph;
	igraph_vector_t iEdges;
	igraph_vector_t iLength;
    
	igraph_vector_init(&iLength,Edges.size());
      igraph_vector_init(&iEdges,2*Edges.size());

	int e = 0;
	double tlength = 0.0;
	for (std::list<_edge>::iterator it = Edges.begin(), end = Edges.end(); it != end; ++it) {
		VECTOR(iEdges)[2*e+0] = (igraph_integer_t) (*it).i;
		VECTOR(iEdges)[2*e+1] = (igraph_integer_t) (*it).j;
		VECTOR(iLength)[e] = (igraph_real_t) (*it).length;
		tlength += (*it).length;
		e++;
	}
	
	igraph_create(&iGraph,&iEdges,N,false);
	SETEANV(&iGraph,"Length",&iLength);

	igraph_vector_destroy(&iEdges);
	igraph_vector_destroy(&iLength);

 	M[0] = igraph_vcount(&iGraph);
	M[1] = igraph_ecount(&iGraph);
	M[2] = tlength;
	M[3] = tlength / M[1];

	int no;
	igraph_vector_t iMem;
	igraph_vector_t iCsize;
	igraph_vector_init(&iMem,0);
	igraph_vector_init(&iCsize,0);

	igraph_clusters(&iGraph,&iMem,&iCsize,&no,IGRAPH_WEAK);

	M[4] = (double)no;

	int smax = 0;
	for (int c = 0; c < no; c++) {
		smax = (VECTOR(iCsize)[c]>smax) ? VECTOR(iCsize)[c] : smax;
	}

	M[5] = (double)smax / M[0];

	M[6] = M[7] = M[8] = M[9] = 0.0;

	int k;
	igraph_vector_t K;
	igraph_vector_init(&K,0);
	igraph_degree(&iGraph,&K,igraph_vss_all(),IGRAPH_ALL,IGRAPH_NO_LOOPS);
	for (int i = 0; i < N; i++) {
		k = (int)VECTOR(K)[i];
		if (k<4) {
			M[5+k] += 1.0/N;
		} else {
			M[9] += 1.0/N;
		}
	}

	igraph_vector_destroy(&iMem);
	igraph_vector_destroy(&iCsize);
	igraph_vector_destroy(&K);

}

void _Graph::GetMinimumDistanceVector(std::vector<double> &Dmin) {
	int i, j;
	double r[3], u[3], d;
	vtkSmartPointer<vtkPoints> Points = vtkSmartPointer<vtkPoints>::New();
	Points -> SetNumberOfPoints(Nodes.size());
	for (i = 0; i < Nodes.size(); i++) {
		Points -> InsertPoint(i,Nodes[i].x,Nodes[i].y,0.0);
	}

	vtkSmartPointer<vtkPolyData> PointsPoly = vtkSmartPointer<vtkPolyData>::New();
	PointsPoly -> SetPoints(Points);

	vtkSmartPointer<vtkKdTreePointLocator> Tree = vtkSmartPointer<vtkKdTreePointLocator>::New();
	Tree -> SetDataSet(PointsPoly);
	Tree -> BuildLocator();

	vtkSmartPointer<vtkIdList> List = vtkSmartPointer<vtkIdList>::New();
	for (i = 0; i < N; i++) {
		Points -> GetPoint(i,r);
		Tree -> FindClosestNPoints(2,r,List);
		Points -> GetPoint(List->GetId(1),u);
		Dmin[i] = sqrt(pow(u[0]-r[0],2)+pow(u[1]-r[1],2));
	}
}

void _Graph::ShuffleCoordinates() {
	for (int i = 0; i < N; i++) {
		Nodes[i].x = Lx * (double)rand() / RAND_MAX;
		Nodes[i].y = Ly * (double)rand() / RAND_MAX;
	}
}

bool _Graph::SwapPairOfEdges(int *io, int *jo, int *po, int *qo) {
 	std::list<_edge>::iterator it1 = Edges.begin();
	std::advance(it1,rand()%Edges.size());
	_edge edge1 = *it1;

	int i = edge1.i;
	int j = edge1.j;
	if (edge1.type < 5) {
		i -= (i/N)*N;
		j -= (j/N)*N;
	} else if (edge1.type == 5) {
		j -= (j/N)*N;
	} else {
		i -= (i/N)*N;
	}

	std::list<_edge>::iterator it2 = Edges.begin();
	std::advance(it2,rand()%Edges.size());
	_edge edge2 = *it2;

	int p = edge2.i;
	int q = edge2.j;
	if (edge2.type < 5) {
		p -= (p/N)*N;
		q -= (q/N)*N;
	} else if (edge2.type == 5) {
		q -= (q/N)*N;
	} else {
		p -= (p/N)*N;
	}

	printf("Edges = (%d,%d)-%d,\t(%d,%d)-%d\n",i,j,edge1.type,p,q,edge2.type);

	ClipNetwork();

	SaveGNET("/Users/matheusviana/GitHub/MItoPlanar/temp.gnet");

	DeleteEdge(edge1);
	//DeleteEdge(edge2);

}

/* =================================================================
   NETWORK MODEL ROUTINES
   =================================================================*/

void GetInstanceOfRandomGraph_Edge(_Graph *Graph, int E, bool _clip) {

	#ifdef DEBUG
		printf("Random Model Constrained by Number of Edges...\n");
		printf("\tAllocating adjacency matrix...\n");
	#endif

	Graph -> CreateVirtualNodes();

	long int k;
	double length;
	int q, i, j, ne = 0;
	int N = Graph -> N;
	vtkSmartPointer<vtkBitArray> ADJ = vtkSmartPointer<vtkBitArray>::New();
	ADJ -> SetNumberOfComponents(1);
	ADJ -> SetNumberOfTuples(N*N);
	ADJ -> FillComponent(0,0);
	for (i=N;i--;) ADJ->SetTuple1(i+i*N,1);

	#ifdef DEBUG
		printf("\tAllocated!\n");
	#endif


	Graph -> Edges.clear();

	while (ne < E) {
		i = rand()%N;
		j = rand()%N;
		if (!ADJ->GetTuple1(i+j*N)) {
			q = Graph -> GetClosestVirtualNode(i,j);
			length = Graph->GetEdgeLength(i,j + q*N);
			if (q) {
				_edge edge = {i,j+q*Graph->N,5,length};
				Graph->Edges.push_back(edge);
				_edge edge_m = {i+_mirror(q)*Graph->N,j,6,length};
				Graph->Edges.push_back(edge_m);
			} else {
				for (q = 0; q < 5; q++) {
		 			_edge edge = {i+q*Graph->N,j+q*Graph->N,q,length};
		 			Graph->Edges.push_back(edge);
		 		}
			}
	 		ne++;
			ADJ -> SetTuple1(i+j*N,1);
			ADJ -> SetTuple1(j+i*N,1);
		}
	}

	if (_clip) Graph -> ClipNetwork();

	#ifdef DEBUG
		printf("\tModel Complete!\n");
	#endif

}

void GetInstanceOfRandomPlanarGraph_Edge(_Graph *Graph, int E, bool _clip) {

	#ifdef DEBUG
		printf("Random Model Constrained by Number of Edges...\n");
		printf("\tAllocating adjacency matrix...\n");
	#endif

	Graph -> CreateVirtualNodes();

	long int k;
	bool connected;
	double length;
	int q, i, j, ne = 0;
	int N = Graph -> N;
	vtkSmartPointer<vtkBitArray> ADJ = vtkSmartPointer<vtkBitArray>::New();
	ADJ -> SetNumberOfComponents(1);
	ADJ -> SetNumberOfTuples(N*N);
	ADJ -> FillComponent(0,0);
	for (i=N;i--;) ADJ->SetTuple1(i+i*N,1);

	#ifdef DEBUG
		printf("\tAllocated!\n");
	#endif


	Graph -> Edges.clear();

	while (ne < E) {
		i = rand()%N;
		j = rand()%N;
		connected = 0;
		if (!ADJ->GetTuple1(i+j*N)) {
			q = Graph -> GetClosestVirtualNode(i,j);
			length = Graph->GetEdgeLength(i,j + q*N);
			if (q) {
				if (Graph->IsPlanarEdge(i,j+q*Graph->N) && Graph->IsPlanarEdge(i+_mirror(q)*Graph->N,j)) {
					_edge edge = {i,j+q*Graph->N,5,length};
					Graph->Edges.push_back(edge);
					_edge edge_m = {i+_mirror(q)*Graph->N,j,6,length};
					Graph->Edges.push_back(edge_m);
					connected = 1;
				}
			} else {
				if (Graph->IsPlanarEdge(i,j)) {
					for (q = 0; q < 5; q++) {
			 			_edge edge = {i+q*Graph->N,j+q*Graph->N,q,length};
			 			Graph->Edges.push_back(edge);
			 		}
			 		connected = 1;
				}
			}
		}
		if (connected) {
	 		ne++;
			ADJ -> SetTuple1(i+j*N,1);
			ADJ -> SetTuple1(j+i*N,1);
		}
	}

	if (_clip) Graph -> ClipNetwork();

	#ifdef DEBUG
		printf("\tModel Complete!\n");
	#endif

}

void GetInstanceOfRandomPlanarGraph_Length(_Graph *Graph, double L, bool _clip) {

	#ifdef DEBUG
		printf("Random Model Constrained by Total Length...\n");
		printf("\tAllocating adjacency matrix...\n");
	#endif

	Graph -> CreateVirtualNodes();

	long int k;
	int q, i, j;
	bool connected;
	int N = Graph -> N;
	double length, total_length = 0.0;
	vtkSmartPointer<vtkBitArray> ADJ = vtkSmartPointer<vtkBitArray>::New();
	ADJ -> SetNumberOfComponents(1);
	ADJ -> SetNumberOfTuples(N*N);
	ADJ -> FillComponent(0,0);
	for (i=N;i--;) ADJ->SetTuple1(i+i*N,1);

	#ifdef DEBUG
		printf("\tAllocated!\n");
	#endif

	Graph -> Edges.clear();

	while (total_length < L) {
		i = rand()%N;
		j = rand()%N;
		connected = 0;
		if (!ADJ->GetTuple1(i+j*N)) {
			q = Graph -> GetClosestVirtualNode(i,j);
			length = Graph->GetEdgeLength(i,j + q*N);
			if (q) {
				if (Graph->IsPlanarEdge(i,j+q*Graph->N) && Graph->IsPlanarEdge(i+_mirror(q)*Graph->N,j)) {
					_edge edge = {i,j+q*Graph->N,5,length};
					Graph->Edges.push_back(edge);
					_edge edge_m = {i+_mirror(q)*Graph->N,j,6,length};
					Graph->Edges.push_back(edge_m);
					connected = 1;
				}
			} else {
				if (Graph->IsPlanarEdge(i,j)) {
					for (q = 0; q < 5; q++) {
			 			_edge edge = {i+q*Graph->N,j+q*Graph->N,q,length};
			 			Graph->Edges.push_back(edge);
			 		}
			 		connected = 1;
				}
			}
		}
		if (connected) {
			ADJ -> SetTuple1(i+j*N,1);
			ADJ -> SetTuple1(j+i*N,1);
			total_length += length;
		}
	}

	if (_clip) Graph -> ClipNetwork();

	#ifdef DEBUG
		printf("\tModel Complete!\n");
	#endif

}

void GetInstanceOfRandomPlanarGraph_Pk1(_Graph *Graph, double pk1, bool _clip) {

	#ifdef DEBUG
		printf("Random Model Constrained by Degree Distribution...\n");
		printf("\tCreating Degree Distribution...\n");
	#endif

	long int k1, k2;
	int q, i, j, pi, pj, temp;
	int N = Graph -> N;
	int nk1 = (int)(pk1*N);
	int nk3 = N - nk1;

	std::vector<int> K (nk1,1);
	K.insert(K.begin(),nk3,3);

	std::vector<int> ID;
	for (i=0;i<N;i++) ID.push_back(i);
	std::random_shuffle(ID.begin(),ID.end(),_rng);

	std::vector<int> V;
	for (i=0;i<N;i++) V.insert(V.begin(),K[i],ID[i]);

	K.clear();
	ID.clear();

	#ifdef DEBUG
		printf("\tDone!\n");
		printf("\tAllocating Adjacency Matrix...\n");
	#endif
	
	Graph -> CreateVirtualNodes();
	
	double length;
	bool _connected;
	vtkSmartPointer<vtkBitArray> ADJ = vtkSmartPointer<vtkBitArray>::New();
	ADJ -> SetNumberOfComponents(1);
	ADJ -> SetNumberOfTuples(N*N);
	ADJ -> FillComponent(0,0);
	for (i=N;i--;) ADJ->SetTuple1(i+i*N,1);

	#ifdef DEBUG
		printf("\tAllocated!\n");
	#endif

	Graph -> Edges.clear();
	Graph -> clock.reset();

	while (V.size() > 1) {
		pi = rand()%V.size();
		pj = rand()%V.size();
		if ( pi > pj ) { temp = pj; pj = pi; pi = temp; }

		i = V[pi];
		j = V[pj];
		_connected = 0;
		if (!ADJ->GetTuple1(i+j*N)) {
			q = Graph -> GetClosestVirtualNode(i,j);
			length = Graph->GetEdgeLength(i,j + q*N);
			if (q) {
				if (Graph->IsPlanarEdge(i,j+q*Graph->N) && Graph->IsPlanarEdge(i+_mirror(q)*Graph->N,j)) {
					_edge edge = {i,j+q*Graph->N,5,length};
					Graph->Edges.push_back(edge);
					_edge edge_m = {i+_mirror(q)*Graph->N,j,6,length};
					Graph->Edges.push_back(edge_m);
					_connected = 1;
				}
			} else {
				if (Graph->IsPlanarEdge(i,j)) {
					for (q = 0; q < 5; q++) {
			 			_edge edge = {i+q*Graph->N,j+q*Graph->N,q,length};
			 			Graph->Edges.push_back(edge);
			 		}
					_connected = 1;
				}
			}
		}
		if (_connected) {
			Graph -> clock.reset();
			ADJ -> SetTuple1(i+j*N,1);
			ADJ -> SetTuple1(j+i*N,1);
			std::swap(V[pj],V.back());
			V.pop_back();
			std::swap(V[pi],V.back());
			V.pop_back();
		} else {
			Graph -> clock.t++;
			if (Graph -> clock.t > N) {
				Graph -> clock.t = 0;
				Graph -> clock.nreset++;
			}
			if (Graph -> clock.nreset > log(N)) {
				Graph -> clock.reset();
				Graph -> DeleteRandomEdge(&i,&j);
				ADJ -> SetTuple1(i+j*N,0);
				ADJ -> SetTuple1(j+i*N,0);
				V.push_back(i);
				V.push_back(j);
			}
		}
	}

	if ( _clip ) Graph -> ClipNetwork();

	#ifdef DEBUG
		printf("\tModel Complete!\n");
	#endif


}

void GetInstanceOfRandomPlanarGraph_Wax(_Graph *Graph, double pk1, double alpha, bool _clip) {

	#ifdef DEBUG
		printf("Random Model Constrained by Degree Distribution...\n");
		printf("\tCreating Degree Distribution...\n");
	#endif

	long int k1, k2;
	int q, i, j, pi, pj, temp;
	int N = Graph -> N;
	int nk1 = (int)(pk1*N);
	int nk3 = N - nk1;

	std::vector<int> K (nk1,1);
	K.insert(K.begin(),nk3,3);

	std::vector<int> ID;
	for (i=0;i<N;i++) ID.push_back(i);
	std::random_shuffle(ID.begin(),ID.end());

	std::vector<int> V;
	for (i=0;i<N;i++) V.insert(V.begin(),K[i],ID[i]);

	K.clear();
	ID.clear();

	#ifdef DEBUG
		printf("\tDone!\n");
		printf("\tAllocating Adjacency Matrix...\n");
	#endif
	
	Graph -> CreateVirtualNodes();
	std::vector<double> Dmin(Graph->N,0.0);
	Graph -> GetMinimumDistanceVector(Dmin);

	bool _connected;
	double dij, pij, pir;
	vtkSmartPointer<vtkBitArray> ADJ = vtkSmartPointer<vtkBitArray>::New();
	ADJ -> SetNumberOfComponents(1);
	ADJ -> SetNumberOfTuples(N*N);
	ADJ -> FillComponent(0,0);
	for (i=N;i--;) ADJ->SetTuple1(i+i*N,1);

	#ifdef DEBUG
		printf("\tAllocated!\n");
	#endif

	Graph -> Edges.clear();
	Graph -> clock.reset();

	while (V.size() > 1) {

		//Accept-reject method
		do {
			pi = rand()%V.size();
			pj = rand()%V.size();
			if ( pi > pj ) { temp = pj; pj = pi; pi = temp; }

			i = V[pi];
			j = V[pj];

			q = Graph -> GetClosestVirtualNode(i,j);
			dij = Graph->GetEdgeLength(i,j + q*N);
			pij = exp(-alpha*dij);
			pir = exp(-alpha*Dmin[i]) * (((double)rand())/RAND_MAX);

		} while (pir>pij);

		_connected = 0;
		if (!ADJ->GetTuple1(i+j*N)) {
			if (q) {
				if (Graph->IsPlanarEdge(i,j+q*Graph->N) && Graph->IsPlanarEdge(i+_mirror(q)*Graph->N,j)) {
					_edge edge = {i,j+q*Graph->N,5,dij};
					Graph->Edges.push_back(edge);
					_edge edge_m = {i+_mirror(q)*Graph->N,j,6,dij};
					Graph->Edges.push_back(edge_m);
					_connected = 1;
				}
			} else {
				if (Graph->IsPlanarEdge(i,j)) {
					for (q = 0; q < 5; q++) {
			 			_edge edge = {i+q*Graph->N,j+q*Graph->N,q,dij};
			 			Graph->Edges.push_back(edge);
			 		}
					_connected = 1;
				}
			}
		}
		if (_connected) {
			Graph -> clock.reset();
			ADJ -> SetTuple1(i+j*N,1);
			ADJ -> SetTuple1(j+i*N,1);
			std::swap(V[pj],V.back()); //Deleting pj first since pj > pi by construction.
			V.pop_back();
			std::swap(V[pi],V.back());
			V.pop_back();
		} else {
			Graph -> clock.t++;
			if (Graph -> clock.t > N) {
				Graph -> clock.t = 0;
				Graph -> clock.nreset++;
			}
			if (Graph -> clock.nreset > log(N)) {
				Graph -> clock.reset();
				Graph -> DeleteRandomEdge(&i,&j);
				ADJ -> SetTuple1(i+j*N,0);
				ADJ -> SetTuple1(j+i*N,0);
				V.push_back(i);
				V.push_back(j);
			}
		}
	}

	if ( _clip ) Graph -> ClipNetwork();

	#ifdef DEBUG
		printf("\tModel Complete!\n");
	#endif


}


/* =================================================================
   MAIN
   =================================================================*/


int main(int argc, char *argv[]) {     

	srand(getpid());

	int _mode = 0;
	int sn, n = 0;
	int nreal = 100;
	bool _save = false;
	double alpha = 1.0;
	double pk1 = 0.75;
	char _Model[4] = {"WAX"};
	char _RootFolder[256] = {""};
	char _SaveFolder[256] = {""};

	for (int i = 0; i < argc; i++) {
		if (!strcmp(argv[i],"-help")) {
			_help();
			return 0;
		}
		if (!strcmp(argv[i],"-path")) {
			sprintf(_RootFolder,"%s//",argv[i+1]);
		}
		if (!strcmp(argv[i],"-model")) {
			sprintf(_Model,"%s",argv[i+1]);
		}
		if (!strcmp(argv[i],"-r")) {
			nreal = atoi(argv[i+1]);
		}
		if (!strcmp(argv[i],"-alpha")) {
			alpha = atof(argv[i+1]);
		}
		if (!strcmp(argv[i],"-pk1")) {
			pk1 = atof(argv[i+1]);
		}
		if (!strcmp(argv[i],"-save")) {
			sprintf(_SaveFolder,"%s",argv[i+1]);
			_save = true;
		}
		if (!strcmp(argv[i],"-pk1_mode")) {
			_mode = 1;
			n = atoi(argv[i+1]);	
			sn = atoi(argv[i+2]);
		}
		if (!strcmp(argv[i],"-size_mode")) {
			_mode = 2;
			n = atoi(argv[i+1]);	
			sn = atoi(argv[i+2]);
		}
		if (!strcmp(argv[i],"-check_planarity")) {
			_mode = 3;
		}
	}

	// Summary file
	char _SUMMFile[256];
	double *M = new double[10];

	sprintf(_SUMMFile,"%s%s.model",_RootFolder,_Model);

	FILE *s = fopen(_SUMMFile,"w");
	fprintf(s,"MitoPlanar V1.0\n");
	if (n) {
		fprintf(s,"Simulation not based on real networks.\n");
		fprintf(s,"Number of nodes: G(%d,%d).\n",n,sn);
	} else {
		fprintf(s,"Folder: %s\n",_RootFolder);
	}
	if (!strcmp(_Model,"WAX")) {
		fprintf(s,"Network Model: %s, alpha = %1.3f\n",_Model,alpha);
	} else {
		fprintf(s,"Network Model: %s\n",_Model);
	}
	fprintf(s,"Number of realizations: %d\n",nreal);
	time_t now = time(0);
	fprintf(s,"%s\n",ctime(&now));
	if (n) {
		fprintf(s,"Pk1\t Model\t NNodes\t NEdges\t Length (um)\t AvgEdgeLength (um)\t NClusters\t Phi\n");
	} else {
		fprintf(s,"Path\t Model\t NNodes\t NEdges\t Length (um)\t AvgEdgeLength (um)\t NClusters\t Phi\n");
	}
	fclose(s);

	if ( _mode==0 ) {

		// Generating list of files to run
		char _cmd[256];
		sprintf(_cmd,"ls %s*.gnet | sed -e 's/.gnet//' > %smitoplanar.files",_RootFolder,_RootFolder);
		system(_cmd);

		int E;
		double L, total_length;

		char _GNETFile[256];
		char _GNETList[256];
		char _MODLList[256];
		char ModelName[256];

		// List of files to run

		sprintf(_GNETList,"%smitoplanar.files",_RootFolder);	
		FILE *f = fopen(_GNETList,"r");

		// Main loop

		while (fgets(_GNETFile,256, f) != NULL) {
			_GNETFile[strcspn(_GNETFile, "\n" )] = '\0';

			printf("%s\n",_GNETFile);

			_Graph Graph;
			Graph.MakeShallowCopy(_GNETFile,&E,&L,&pk1);

			total_length = 0.0;

			for (int net = 0; net <  nreal; net++) {

				if (!strcmp(_Model,"EDG")) {
					GetInstanceOfRandomPlanarGraph_Edge(&Graph,E,true);
				}
				if (!strcmp(_Model,"LGT")) {
					GetInstanceOfRandomPlanarGraph_Length(&Graph,L,true);
				}
				if (!strcmp(_Model,"PK1")) {
					GetInstanceOfRandomPlanarGraph_Pk1(&Graph,pk1,true);
				}
				if (!strcmp(_Model,"WAX")) {
					GetInstanceOfRandomPlanarGraph_Wax(&Graph,pk1,alpha,true);
				}
				if (!strcmp(_Model,"WXS")) {
					Graph.ShuffleCoordinates();
					GetInstanceOfRandomPlanarGraph_Wax(&Graph,pk1,alpha,true);
				}
				Graph.GetProperties(M);

				s = fopen(_SUMMFile,"a");
				fprintf(s,"%s\t%s\t%d\t%d\t%1.3f\t%1.3f\t%d\t%1.3f\n",_GNETFile,_Model,(int)M[0],(int)M[1],M[2],M[3],(int)M[4],M[5]);
				fclose(s);

				total_length += M[2];

				sprintf(_MODLList,"%s-%d.%s",_GNETFile,net,_Model);
				if (_save) Graph.SaveGNET(_MODLList);

				printf("\t%d\n",net);
			}

			printf("\tTotal length = %1.5f\n",total_length/nreal);
		
		}

		fclose(f);

	} else if (_mode == 1) {

		#ifdef DEBUG
			printf("Simulation mode (pk1 mode)...\n");
		#endif

		double pk1 = 0.05;

		do {

			for (int net = 0; net <  nreal; net++) {

				_Graph Graph;
				Graph.GenerateFromScratch(n,sn);

				if (!strcmp(_Model,"WAX")) {
					GetInstanceOfRandomPlanarGraph_Wax(&Graph,pk1,alpha,true);
				}

				Graph.GetProperties(M);

				s = fopen(_SUMMFile,"a");
				fprintf(s,"%1.2f\t%s\t%d\t%d\t%1.3f\t%1.3f\t%d\t%1.3f\n",pk1,_Model,(int)M[0],(int)M[1],M[2],M[3],(int)M[4],M[5]);
				fclose(s);

			}

			pk1 += 0.05;

		} while (pk1 < 1.01);

	} else if (_mode==2) {

		#ifdef DEBUG
			printf("Simulation mode (size mode)...\n");
		#endif

		for (int nnodes = n; nnodes <= 300; nnodes += sn) {

			for (int net = 0; net < nreal; net++) {

				_Graph Graph;
				Graph.GenerateFromScratch(nnodes,0);

				if (!strcmp(_Model,"WAX")) {
					GetInstanceOfRandomPlanarGraph_Wax(&Graph,pk1,alpha,true);
				}

				Graph.GetProperties(M);

				s = fopen(_SUMMFile,"a");
				fprintf(s,"0.75\t%s\t%d\t%d\t%1.3f\t%1.3f\t%d\t%1.3f\n",_Model,(int)M[0],(int)M[1],M[2],M[3],(int)M[4],M[5]);
				fclose(s);

			}

			printf("N = %d\n",nnodes);

		}

	} else if (_mode==3) {

		#ifdef DEBUG
			printf("Checking for planarity...\n");
		#endif

		// Generating list of files to run
		char _cmd[256];
		sprintf(_cmd,"ls %s*.gnet | sed -e 's/.gnet//' > %smitoplanar.files",_RootFolder,_RootFolder);
		system(_cmd);

		int E, ncr;
		double L, pk1, total_length;

		char _GNETFile[256];
		char _GNETList[256];
		char _RANDName[256];

		// List of files to run

		sprintf(_GNETList,"%smitoplanar.files",_RootFolder);	
		FILE *f = fopen(_GNETList,"r");

		while (fgets(_GNETFile,256, f) != NULL) {
			_GNETFile[strcspn(_GNETFile, "\n" )] = '\0';

			printf("%s\n",_GNETFile);

			sprintf(_RANDName,"%s.random",_GNETFile);
			FILE *fs = fopen(_RANDName,"w");
			fprintf(fs,"Type\tN\tE\tL\t<l>\tNc\tPhi\tNcr\tP1\tP2\tP3\tP4+\n");

			_Graph Graph;
			Graph.MakeDeepCopy(_GNETFile,&E,&L,&pk1);

			ncr = Graph.GetTotalNumberOfCrosses();
			Graph.ClipNetwork();
			Graph.GetProperties(M);
			fprintf(fs,"Real\t%d\t%d\t%1.3f\t%1.3f\t%d\t%1.3f\t%d\t%1.2f\t%1.2f\t%1.2f\t%1.2f\n",(int)M[0],(int)M[1],M[2],M[3],(int)M[4],M[5],ncr,M[6],M[7],M[8],M[9]);

			for (int net = 0; net <  nreal; net++) {

				GetInstanceOfRandomGraph_Edge(&Graph,E,false);

				ncr = Graph.GetTotalNumberOfCrosses();
				Graph.ClipNetwork();
				Graph.GetProperties(M);
				fprintf(fs,"Rand\t%d\t%d\t%1.3f\t%1.3f\t%d\t%1.3f\t%d\t%1.2f\t%1.2f\t%1.2f\t%1.2f\n",(int)M[0],(int)M[1],M[2],M[3],(int)M[4],M[5],ncr,M[6],M[7],M[8],M[9]);

			}

			fclose(fs);

		}

		fclose(f);

	}

	return 0;
}