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

struct _node{
	double x;
	double y;
};

struct _edge{
	int i;
	int j;
	int type;
	double length;
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
	double xmin, ymin, xmax, ymax;
	std::vector<_node> Nodes;
	std::list<_edge> Edges;
	_Graph() {
		N = 0;
	};
	void CreateVirtualNodes();
	bool IsPlanarEdge(int i, int j);
	void SavePolyData(const char FileName[]);
	int GetClosestVirtualNode(int i, int j);
	double GetEdgeLength(int i, int j);
	void MakeShallowCopy(const char FilePrefix[], int *NEdges, double *Length);
	void ClipNetwork();
	void SaveGNET(const char FileName[]);
};

void _Graph::CreateVirtualNodes() {
	double r[3], x, y;
	N = Nodes.size();
	double Q[4][2] = {{0,ymax},{xmax,0},{0,-ymax},{-xmax,0}};
	for (int q = 1; q <= 4; q++) {
		for (int i = 0; i < N; i++) {
			_node node = {Nodes[i].x+Q[q-1][0],Nodes[i].y+Q[q-1][1]};
			Nodes.push_back(node);
		}
	}
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
	    	Q =     ri[0] * ( rv[1] - ru[1] ) + ru[0] * ( ri[1] - rv[1] ) + rv[0] * ( ru[1] - ri[1] );
			xc = ri[0] + (Q / Det) * ( rj[0] - ri[0] );
			yc = ri[1] + (Q / Det) * ( rj[1] - ri[1] );

			if (SegmentsContainPoint(ri,rj,ru,rv,xc,yc)) return false;
        }
    }
	return true;
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

double _Graph::GetEdgeLength(int i, int j) {
	double ri[3], rj[3];
	ri[0] = Nodes[i].x;
	ri[1] = Nodes[i].y;
	rj[0] = Nodes[j].x;
	rj[1] = Nodes[j].y;
	return sqrt(pow(ri[0]-rj[0],2)+pow(ri[1]-rj[1],2));
}

void _Graph::MakeShallowCopy(const char FilePrefix[], int *NEdges, double *Length) {

	#ifdef DEBUG
		printf("Copying Coordinates From %s.coo\n",FilePrefix);
	#endif

    char _path[128];
	int i, j, E = 0;
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

	Nodes.clear();

	#ifdef DEBUG
		printf("\t#Nodes = %d\n",N);
		printf("\t#Edges = %d\n",E);
		printf("\t#Length = %1.3f\n",L);
	#endif

	xmax = ymax = 0.0;
	xmin = ymin = 1E6;
	while (fscanf(fc,"%f %f %f",&x,&y,&z)!=EOF) {
		x += _eps * (1.0*rand())/RAND_MAX;
		y += _eps * (1.0*rand())/RAND_MAX;
		_node node = {x,y};
		Nodes.push_back(node);
		xmin = (x<xmin) ? x : xmin;
		xmax = (x>xmax) ? x : xmax;
		ymin = (y<ymin) ? y : ymin;
		ymax = (y>xmax) ? y : xmax;
	}
	fclose(fc);
	*Length = L;
	*NEdges = E;

	#ifdef DEBUG
		printf("\t%1.3f<=x<=%1.3f\n",xmin,xmax);
		printf("\t%1.3f<=y<=%1.3f\n",ymin,ymax);
		printf("\tCopy Complete!\n");
	#endif
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


	Graph -> Edges.clear();

	while (ne < E) {
		i = rand()%N;
		j = rand()%N;
		k = i + j *N;
		if (!ADJ->GetTuple1(k)) {
			q = Graph -> GetClosestVirtualNode(i,j);
			if (q) {
				if (Graph->IsPlanarEdge(i,j+q*Graph->N) && Graph->IsPlanarEdge(i+_mirror(q)*Graph->N,j)) {
					_edge edge = {i,j+q*Graph->N,5,0.0};
					Graph->Edges.push_back(edge);
					_edge edge_m = {i+_mirror(q)*Graph->N,j,6,0.0};
					Graph->Edges.push_back(edge_m);
					ne++;
					ADJ -> SetTuple1(k,1);
					ADJ -> SetTuple1(k,1);
				}
			} else {
				if (Graph->IsPlanarEdge(i,j)) {
					for (q = 0; q < 5; q++) {
			 			_edge edge = {i+q*Graph->N,j+q*Graph->N,q,0.0};
			 			Graph->Edges.push_back(edge);
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

	Graph -> Edges.clear();

	while (total_length < L) {
		i = rand()%N;
		j = rand()%N;
		k = i + j *N;
		if (!ADJ->GetTuple1(k)) {
			q = Graph -> GetClosestVirtualNode(i,j);
			if (q) {
				if (Graph->IsPlanarEdge(i,j+q*Graph->N) && Graph->IsPlanarEdge(i+_mirror(q)*Graph->N,j)) {
					_edge edge = {i,j+q*Graph->N,5,0.0};
					Graph->Edges.push_back(edge);
					_edge edge_m = {i+_mirror(q)*Graph->N,j,6,0.0};
					Graph->Edges.push_back(edge_m);
					ADJ -> SetTuple1(k,1);
					ADJ -> SetTuple1(k,1);
					total_length += Graph -> GetEdgeLength(i,j+q*Graph->N);
				}
			} else {
				if (Graph->IsPlanarEdge(i,j)) {
					for (q = 0; q < 5; q++) {
			 			_edge edge = {i+q*Graph->N,j+q*Graph->N,q,0.0};
			 			Graph->Edges.push_back(edge);
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

void _Graph::ClipNetwork() {
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
	for (i=N;i<5*N;i++) Nodes.erase(Nodes.end());
}

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

int main() {

	srand(getpid());

	int E;
	double L;
	char OriginalName[128], ModelName[128];
	sprintf(OriginalName,"101010_1");

	_Graph Graph;
	Graph.MakeShallowCopy(OriginalName,&E,&L);

	for (int net = 0; net <  5; net++) {

		Graph.CreateVirtualNodes();
		GetInstanceOfRandomPlanarGraph_Edge(&Graph,E);
		Graph.ClipNetwork();

		sprintf(ModelName,"%s-%03d.vtk",OriginalName,net);
		Graph.SavePolyData(ModelName);

		sprintf(ModelName,"%s-%03d.gnet",OriginalName,net);
		Graph.SaveGNET(ModelName);

	}
	return 0;
}