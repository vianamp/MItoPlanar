#include <list>
#include <omp.h>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <igraph/igraph.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkCellArray.h>
#include <vtkBitArray.h>

#define _min(a,b) ((a<b)?a:b)
#define _max(a,b) ((a>b)?a:b)
#define _eps 1E-6

struct _edge{
	int i;
	int j;
	bool frontier;
};

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

class _Graph{
private:
public:
	int N;
	vtkSmartPointer<vtkPoints> Nodes;
	std::list<_edge> LEdges;
	_Graph() {
		Nodes = vtkSmartPointer<vtkPoints>::New();		
	}
	void CreateVirtualNodes();
	bool IsPlanarEdge(int i, int j);
	void SavePolyData(const char FileName[]);
	int GetClosestVirtualNode(int i, int j);
	double GetEdgeLength(int i, int j);
};

void _Graph::CreateVirtualNodes() {
	double r[3];
	N = Nodes -> GetNumberOfPoints();
	double *B = Nodes -> GetBounds();
	double Q[4][2] = {{0,B[3]},{B[1],0},{0,-B[3]},{-B[1],0}};
	for (int q = 1; q <= 4; q++) {
		for (int i = 0; i < N; i++) {
			Nodes -> GetPoint(i,r);
			Nodes -> InsertPoint(i+q*N,r[0]+Q[q-1][0],r[1]+Q[q-1][1],0.0);
		}
	}
}

bool _Graph::IsPlanarEdge(int i, int j) {
	double Det, Q, xc, yc, fac;
	double ri[3], rj[3], ru[3], rv[3];
	Nodes -> GetPoint(i,ri);
	Nodes -> GetPoint(j,rj);

	for (std::list<_edge>::const_iterator it = LEdges.begin(), end = LEdges.end(); it != end; ++it) {
		Nodes -> GetPoint((*it).i,ru);
		Nodes -> GetPoint((*it).j,rv);

		Det = ( rj[1] - ri[1] ) * ( rv[0] - ru[0] ) - ( rv[1] - ru[1] ) * ( rj[0] - ri[0] );

	    if(fabs(Det) > _eps) {
	    	Q =     ri[0] * ( rv[1] - ru[1] ) + ru[0] * ( ri[1] - rv[1] ) + rv[0] * ( ru[1] - ri[1] );
			xc = ri[0] + (Q / Det) * ( rj[0] - ri[0] );
			yc = ri[1] + (Q / Det) * ( rj[1] - ri[1] );

			if (SegmentsContainPoint(ri,rj,ru,rv,xc,yc)) return false;
        }
    }
	return true;
}

void _Graph::SavePolyData(const char FileName[]) {
	printf("Saving polydata...\n");
	vtkSmartPointer<vtkCellArray> Cells = vtkSmartPointer<vtkCellArray>::New();
	for (std::list<_edge>::const_iterator it = LEdges.begin(), end = LEdges.end(); it != end; ++it) {
		Cells -> InsertNextCell(2);
		Cells -> InsertCellPoint((*it).i);
		Cells -> InsertCellPoint((*it).j);
	}

	vtkSmartPointer<vtkPolyData> Net = vtkSmartPointer<vtkPolyData>::New();
	Net -> SetPoints(Nodes);
	Net -> SetLines(Cells);

	vtkSmartPointer<vtkPolyDataWriter> PolyWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
	PolyWriter -> SetFileName(FileName);
	PolyWriter -> SetInputData(Net);
	PolyWriter -> Write();
}

int _Graph::GetClosestVirtualNode(int i, int j) {
	int qo;
	double ri[3], rj[3], d, dmin = 1E5;
	Nodes -> GetPoint(i,ri);
	for (int q = 0; q < 5; q++) {
		Nodes -> GetPoint(j+q*N,rj);
		d = pow(rj[1]-ri[1],2) + pow(rj[0]-ri[0],2);
		if (d < dmin) {
			dmin = d;
			qo = q;
		}
	}
	return qo;
}

double _Graph::GetEdgeLength(int i, int j) {
	double ri[3], rj[3];
	Nodes -> GetPoint(i,ri);
	Nodes -> GetPoint(j,rj);
	return sqrt(pow(ri[0]-rj[0],2)+pow(ri[1]-rj[1],2));
}

void GetInstanceOfRandomPlanarGraph_Edge(_Graph *Graph, int E) {

	#ifdef DEBUG
		printf("Random Model Constrained by Total Length...\n");
		printf("\tAllocating adjacency matrix...\n");
	#endif

	long int k;
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


	Graph -> LEdges.clear();

	while (ne < E) {
		i = rand()%N;
		j = rand()%N;
		k = i + j *N;
		if (!ADJ->GetTuple1(k)) {
			q = Graph -> GetClosestVirtualNode(i,j);
			if (q) {
				if (Graph->IsPlanarEdge(i,j+q*Graph->N) && Graph->IsPlanarEdge(i+_mirror(q)*Graph->N,j)) {
					_edge edge = {i,j+q*Graph->N,1};
					Graph->LEdges.push_back(edge);
					_edge edge_m = {i+_mirror(q)*Graph->N,j,1};
					Graph->LEdges.push_back(edge_m);
					ne++;
					ADJ -> SetTuple1(k,1);
					ADJ -> SetTuple1(k,1);
				}
			} else {
				if (Graph->IsPlanarEdge(i,j)) {
					for (q = 0; q < 5; q++) {
			 			_edge edge = {i+q*Graph->N,j+q*Graph->N,0};
			 			Graph->LEdges.push_back(edge);
			 		}
			 		ne++;
					ADJ -> SetTuple1(k,1);
					ADJ -> SetTuple1(k,1);
				}
			}
		}
	}

	#ifdef DEBUG
		printf("\tModel Complete!\n");
	#endif

}

void GetInstanceOfRandomPlanarGraph_Length(_Graph *Graph, double L) {

	#ifdef DEBUG
		printf("Random Model Constrained by Total Length...\n");
		printf("\tAllocating adjacency matrix...\n");
	#endif

	long int k;
	int q, i, j;
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

	Graph -> LEdges.clear();

	while (total_length < L) {
		i = rand()%N;
		j = rand()%N;
		k = i + j *N;
		if (!ADJ->GetTuple1(k)) {
			q = Graph -> GetClosestVirtualNode(i,j);
			if (q) {
				if (Graph->IsPlanarEdge(i,j+q*Graph->N) && Graph->IsPlanarEdge(i+_mirror(q)*Graph->N,j)) {
					_edge edge = {i,j+q*Graph->N,1};
					Graph->LEdges.push_back(edge);
					_edge edge_m = {i+_mirror(q)*Graph->N,j,1};
					Graph->LEdges.push_back(edge_m);
					ADJ -> SetTuple1(k,1);
					ADJ -> SetTuple1(k,1);
					total_length += Graph -> GetEdgeLength(i,j+q*Graph->N);
				}
			} else {
				if (Graph->IsPlanarEdge(i,j)) {
					for (q = 0; q < 5; q++) {
			 			_edge edge = {i+q*Graph->N,j+q*Graph->N,0};
			 			Graph->LEdges.push_back(edge);
			 		}
					ADJ -> SetTuple1(k,1);
					ADJ -> SetTuple1(k,1);
					total_length += Graph -> GetEdgeLength(i,j);
				}
			}
		}
	}

	#ifdef DEBUG
		printf("\tModel Complete!\n");
	#endif

}

_Graph LoadShallowCopy(const char FilePrefix[], int *NEdges, double *Length) {

	#ifdef DEBUG
		printf("Copying Coordinates From %s.coo\n",FilePrefix);
	#endif

	_Graph Graph;
    char _path[128];
	int N, i, j, E = 0;
	float x, y, z, length, L = 0.0;
    sprintf(_path,"%s.gnet",FilePrefix);
	FILE *fg = fopen(_path,"r");
    sprintf(_path,"%s.coo",FilePrefix);
	FILE *fc = fopen(_path,"r");
	fscanf(fg,"%d",&N);
	while (fscanf(fg,"%d %d %f",&i,&j,&length)!=EOF) {
		E++;
		L+=length;
	}
	fclose(fg);

	#ifdef DEBUG
		printf("\t#Nodes = %d\n",N);
		printf("\t#Edges = %d\n",E);
		printf("\t#Length = %1.3f\n",L);
	#endif

	while (fscanf(fc,"%f %f %f",&x,&y,&z)!=EOF) {
		x += _eps * (1.0*rand())/RAND_MAX;
		y += _eps * (1.0*rand())/RAND_MAX;
		Graph.Nodes -> InsertNextPoint(x,y,z);
	}
	fclose(fc);
	*Length = L;
	*NEdges = E;
	Graph.CreateVirtualNodes();

	#ifdef DEBUG
		printf("\tCopy Complete!\n");
	#endif

	return Graph;
}

int main() {

	srand(getpid());

	int E;
	double L;

	_Graph Graph = LoadShallowCopy("101010_1",&E,&L);

	GetInstanceOfRandomPlanarGraph_Edge(&Graph,E);
	Graph.SavePolyData("Net1.vtk");

	GetInstanceOfRandomPlanarGraph_Length(&Graph,20);
	Graph.SavePolyData("Net2.vtk");

	return 0;
}