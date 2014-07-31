#include <list>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkCellArray.h>
#include <vtkBitArray.h>

#define _min(a,b) ((a<b)?a:b)
#define _max(a,b) ((a>b)?a:b)
#define _eps 1E-6

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
	void MakeShallowCopy(const char FilePrefix[], int *NEdges, double *Length);

	double GetEdgeLength(int i, int j);
	bool IsPlanarEdge(int i, int j);
	void CreateVirtualNodes();
	int GetClosestVirtualNode(int i, int j);
	void ClipNetwork();
	void DeleteRandomEdge(int *, int *j);
	
};

/* =================================================================
   AUXILIAR FUNCTIONS
   =================================================================*/

int _mirror(int q) {
	switch (q) {
		case 1: return 3;
		case 2: return 4;
		case 3: return 1;
	 default  : return 2;
	}
}

bool IsInBetween(double ri[3], double rj[3], double v, int d) {
	return ( v > _min(ri[d],rj[d]) ) && ( v < _max(ri[d],rj[d]) ) ? true : false;
}

bool SegmentsContainPoint(double ri[3], double rj[3], double ru[3], double rv[3], double xc, double yc) {
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

/* =================================================================
   I/O METHODS
   =================================================================*/

void _Graph::SavePolyData(const char FileName[]) {
	printf("Saving polydata...\n");

	vtkSmartPointer<vtkPoints> Points = vtkSmartPointer<vtkPoints>::New();
	for (int i=0; i < 5*N; i++) {
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

void _Graph::MakeShallowCopy(const char FilePrefix[], int *NEdges, double *Length) {

	#ifdef DEBUG
		printf("Copying Coordinates From %s.coo2d\n",FilePrefix);
	#endif

	char text[128], _path[128];
	int i, j, E = 0;
	float x, y, length, L = 0.0;
	sprintf(_path,"%s.gnet",FilePrefix);
	FILE *fg = fopen(_path,"r");
	sprintf(_path,"%s.coo2d",FilePrefix);
	FILE *fc = fopen(_path,"r");
	fscanf(fg,"%d",&N);
	while (fscanf(fg,"%d %d %f",&i,&j,&length)!=EOF) {
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


	#ifdef DEBUG
		printf("\tLx = %1.3f, Ly = %1.3f\n",Lx,Ly);
		printf("\tCopy Complete!\n");
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
	double Det, Q, xc, yc, fac;
	double ri[3], rj[3], ru[3], rv[3];
	ri[0] = Nodes[i].x; ri[1] = Nodes[i].y;
	rj[0] = Nodes[j].x; rj[1] = Nodes[j].y;

	for (std::list<_edge>::iterator it = Edges.begin(), end = Edges.end(); it != end; ++it) {
		ru[0] = Nodes[(*it).i].x; ru[1] = Nodes[(*it).i].y;
		rv[0] = Nodes[(*it).j].x; rv[1] = Nodes[(*it).j].y;

		Det = ( rj[1] - ri[1] ) * ( rv[0] - ru[0] ) - ( rv[1] - ru[1] ) * ( rj[0] - ri[0] );

		if(fabs(Det) > _eps) {
			Q =  ri[0] * ( rv[1] - ru[1] ) + ru[0] * ( ri[1] - rv[1] ) + rv[0] * ( ru[1] - ri[1] );
			xc = ri[0] + (Q / Det) * ( rj[0] - ri[0] );
			yc = ri[1] + (Q / Det) * ( rj[1] - ri[1] );

			if (SegmentsContainPoint(ri,rj,ru,rv,xc,yc)) return false;
		}
	}
	return true;
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

void _Graph::DeleteRandomEdge(int *io, int *jo) {

	#ifdef DEBUG
		printf("Deleteting edge...\n");
	#endif

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

/* =================================================================
   NETWORK MODEL ROUTINES
   =================================================================*/

void GetInstanceOfRandomPlanarGraph_Edge(_Graph *Graph, int E) {

	#ifdef DEBUG
		printf("Random Model Constrained by Number of Edges...\n");
		printf("\tAllocating adjacency matrix...\n");
	#endif

	Graph -> CreateVirtualNodes();

	long int k;
	bool connected;
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
			if (q) {
				if (Graph->IsPlanarEdge(i,j+q*Graph->N) && Graph->IsPlanarEdge(i+_mirror(q)*Graph->N,j)) {
					_edge edge = {i,j+q*Graph->N,5,0.0};
					Graph->Edges.push_back(edge);
					_edge edge_m = {i+_mirror(q)*Graph->N,j,6,0.0};
					Graph->Edges.push_back(edge_m);
					connected = 1;
				}
			} else {
				if (Graph->IsPlanarEdge(i,j)) {
					for (q = 0; q < 5; q++) {
			 			_edge edge = {i+q*Graph->N,j+q*Graph->N,q,0.0};
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

	Graph -> ClipNetwork();

	#ifdef DEBUG
		printf("\tModel Complete!\n");
	#endif

}

void GetInstanceOfRandomPlanarGraph_Length(_Graph *Graph, double L) {

	#ifdef DEBUG
		printf("Random Model Constrained by Total Length...\n");
		printf("\tAllocating adjacency matrix...\n");
	#endif

	Graph -> CreateVirtualNodes();

	long int k;
	int q, i, j;
	bool connected;
	int N = Graph -> N;
	double total_length;
	vtkSmartPointer<vtkBitArray> ADJ = vtkSmartPointer<vtkBitArray>::New();
	ADJ -> SetNumberOfComponents(1);
	ADJ -> SetNumberOfTuples(N*N);
	ADJ -> FillComponent(0,0);
	for (i=N;i--;) ADJ->SetTuple1(i+i*N,1);

	#ifdef DEBUG
		printf("\tAllocated!\n");
	#endif

	Graph -> Edges.clear();

	connected = 0;
	while (total_length < L) {
		i = rand()%N;
		j = rand()%N;
		if (!ADJ->GetTuple1(i+j*N)) {
			q = Graph -> GetClosestVirtualNode(i,j);
			if (q) {
				if (Graph->IsPlanarEdge(i,j+q*Graph->N) && Graph->IsPlanarEdge(i+_mirror(q)*Graph->N,j)) {
					_edge edge = {i,j+q*Graph->N,5,0.0};
					Graph->Edges.push_back(edge);
					_edge edge_m = {i+_mirror(q)*Graph->N,j,6,0.0};
					Graph->Edges.push_back(edge_m);
					connected = 1;
				}
			} else {
				if (Graph->IsPlanarEdge(i,j)) {
					for (q = 0; q < 5; q++) {
			 			_edge edge = {i+q*Graph->N,j+q*Graph->N,q,0.0};
			 			Graph->Edges.push_back(edge);
			 		}
			 		connected = 1;
				}
			}
		}
		if (connected) {
			ADJ -> SetTuple1(i+j*N,1);
			ADJ -> SetTuple1(j+i*N,1);
			total_length += Graph -> GetEdgeLength(i,j);			
		}
	}

	Graph -> ClipNetwork();

	#ifdef DEBUG
		printf("\tModel Complete!\n");
	#endif

}

void GetInstanceOfRandomPlanarGraph_Pk1(_Graph *Graph, double pk1) {

	#ifdef DEBUG
		printf("Random Model Constrained by Degree Distribution...\n");
		printf("\tCreating Degree Distribution...\n");
	#endif

	long int k1, k2;
	int q, i, j, pi, pj;
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
	
	bool _connected;
	double total_length;
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
		i = V[pi];
		j = V[pj];
		_connected = 0;
		if (!ADJ->GetTuple1(i+j*N)) {
			q = Graph -> GetClosestVirtualNode(i,j);
			if (q) {
				if (Graph->IsPlanarEdge(i,j+q*Graph->N) && Graph->IsPlanarEdge(i+_mirror(q)*Graph->N,j)) {
					_edge edge = {i,j+q*Graph->N,5,0.0};
					Graph->Edges.push_back(edge);
					_edge edge_m = {i+_mirror(q)*Graph->N,j,6,0.0};
					Graph->Edges.push_back(edge_m);
					_connected = 1;
				}
			} else {
				if (Graph->IsPlanarEdge(i,j)) {
					for (q = 0; q < 5; q++) {
			 			_edge edge = {i+q*Graph->N,j+q*Graph->N,q,0.0};
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
			total_length += Graph -> GetEdgeLength(i,j);
			std::swap(V[pi],V.back());
			V.pop_back();
			std::swap(V[pj],V.back());
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

	Graph -> ClipNetwork();

	#ifdef DEBUG
		printf("\tModel Complete!\n");
	#endif


}


/* =================================================================
   MAIN
   =================================================================*/


int main(int argc, char *argv[]) {     

	srand(getpid());
	char _RootFolder[256] = {""};
	//sprintf(_RootFolder,"");

	for (int i = 0; i < argc; i++) {
		if (!strcmp(argv[i],"-path")) {
			sprintf(_RootFolder,"%s//",argv[i+1]);
		}
	}

	// Generating list of files to run
	char _cmd[256];
	sprintf(_cmd,"ls %s*.gnet | sed -e 's/.gnet//' > %smitoplanar.files",_RootFolder,_RootFolder);
	system(_cmd);

	int E;
	double L;

	char _GNETFile[256];
	char _GNETList[256];
	char ModelName[256];
	sprintf(_GNETList,"%smitoplanar.files",_RootFolder);
	
	FILE *f = fopen(_GNETList,"r");
	while (fgets(_GNETFile,256, f) != NULL) {
		_GNETFile[strcspn(_GNETFile, "\n" )] = '\0';

		_Graph Graph;
		Graph.MakeShallowCopy(_GNETFile,&E,&L);

		for (int net = 0; net <  100; net++) {

			printf("%d\n",net);

			GetInstanceOfRandomPlanarGraph_Edge(&Graph,20);

			//GetInstanceOfRandomPlanarGraph_Length(&Graph,L);

			//GetInstanceOfRandomPlanarGraph_Pk1(&Graph,0.4);

			sprintf(ModelName,"%s-EDG-%03d.vtk",_GNETFile,net);
			Graph.SavePolyData(ModelName);

			//sprintf(ModelName,"%s-EDG-%03d.gnet",_GNETFile,net);
			//Graph.SaveGNET(ModelName);

		}

	}
	return 0;
}