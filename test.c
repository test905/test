
//gcc -Wpedantic -Wall -Wextra -Wno-unused -Wno-unused-parameter -Wno-unused-result -Werror -std=c99 -march=core2 -O3 -pthread cluster.c -o cluster && time ./cluster


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <pthread.h>
#include <fcntl.h>
#include <unistd.h> error


//  define constants


#define B  33554432             // max. file read size + 1
#define N  65536                // max. number of points + 1
#define D  16                   // max. dimension count
#define M  64                   // max. partition count
#define T  8                    // number of threads
#define E  4096                 // max. number of edges + 1
#define P  8                    // decimal precision


//  define flags


//#define PREFETCH


//  define types


typedef  unsigned char   u8;
typedef  unsigned short  u16;
typedef  float           f32;


//  define macros


                                // align to cache-line or page
#define align(x)  __attribute__ ((aligned(x)))
#define a128  align(128)
#define a4m  align(4194304)


#define z1  = {0}               // initialize to zero
#define z3  = {{{0}}}           // 


//  data segment  --------------------------------------------------------------


//  parameter


int variant;

int delta_m;
int tau_factor;
int eps_factor;

float delta;
float tau;

char *file1;
char *file2;


//  initialize


a4m  u8  buffer[B] z1;          // (32 m)  file read buffer

a128 int point_count = N - 1;   //         number of points
a4m  f32 point[N][D];           // (4 m)   data from file
a4m  u16 point_node[N];         // (128 k) point -> node


//  thread control


a4m  u16 task[T + 1];           //         partition work load

a128 pthread_t thread[T];


//  parallel read


a128 int dim_count = D;         //         dimension count

a128 int node_count = 1;        //         number of nodes + 1
a4m  u8  coord[N][D];           // (1 m)   node coords
a4m  u16 density_scaled[N] z1;  // (128 k) number of points in node
a4m  f32 estimate[N][D][2] z3;  // (8 m)   estimate points in node

a4m  u16 index[D][N][M] z3;     // (128 m) coord -> node


//  parallel write


a4m  u16 edge_count[N] z1;      // (128 k) number of edges
a4m  u16 edge[N][E];            // (512 m) edges to adjacent nodes


//  sequential


a128 int rho_scaled;

a128 int stack_count = 0;       //         stack size
a4m  u16 stack[N];              // (128 k) depth-first search

a4m  u16 node_component[N];     // (128 k) node -> component

a128 int component_count;       //         number of components + 1
a4m  u16 component_cluster[N];  // (128 k) component -> cluster

a128 int cluster_count = 1;     //         number of clusters


//  stage 1: read  -------------------------------------------------------------


int clamp(int min,int max,int a) {
	
	
	// bound a to min...max
	//
	//
	// return min if a < min,
	// return max if a > max,
	// return a otherwise
	
	
	if(a < min) {
		a = min;
	}
	else if(a > max) {
		a = max;
	}
	
	return a;
}


int skip(char** cur,char expected) {
	
	
	// increment cursor if character
	// found, otherwise do nothing
	//
	//
	// return 1 if character found,
	// return 0 otherwise
	
	
	if(**cur == expected) {
		++(*cur);
		return 1;
	}
	
	return 0;
}


int convert(char** cur,float* value) {
	
	
	// convert decimal literal to ieee 754
	// format: "9.999"... with optional "e-99"
	//
	//
	// expecting: '0'...'9', 'e', '-', '.'
	// abort on unexpected character
	//
	//
	// return 0 if first character unexpected,
	// return 1 otherwise
	
	
	float exp = 10.f;
	*value = 0.f;
	
	int k = 0;
	for(; k < P; ++k) {
		
		char ch = **cur;
		if(ch == 'e') {
			
			*value = 0.f;
			break;
		}
		else if(ch == '-') {
		}
		else if(ch != '.') {
			
			u8 digit = ch - '0';
			if(digit > 9) {
				break;
			}
			
			exp *= 0.1f;
			*value += exp * digit;
		}
		
		++(*cur);
	}
	for(; k < 32; ++k) {
		
		char ch = **cur;
		if(ch == 'e') {
			
			*value = 0.f;
		}
		else if(ch == '-') {
		}
		else if(ch != '.') {
			
			u8 digit = ch - '0';
			if(digit > 9) {
				break;
			}
		}
		
		++(*cur);
	}
	
	return !!k;
}


void stage1_read() {
	
	
	// stage 1: read input file and initialize data segment
	//
	//
	// read points from file containing comma-separated values
	// identify points with specified partition of input domain
	//
	//
	// count points in node and estimate contained points
	//
	//
	// create index structure to determine node at coordinates
	
	
	int fd1 = open(file1,O_RDONLY);
	
	if(fd1 == -1) {
		fprintf(stderr,"Specified file does not exist or failed to open.\n");
		exit(1);
	}
	
	int size = read(fd1,buffer,sizeof(buffer));
	buffer[sizeof(buffer) - 1] = 0;
	
	close(fd1);
	
	
	char *cur = (char*) buffer;
	char *max = (char*) (buffer + size);
	
	for(int i = 0; i < point_count; ++i) {
		
		for(int j = 0; j < dim_count; ++j) {
			
			skip(&cur,',');
			skip(&cur,' ');
			
			if(!convert(&cur,&point[i][j])) {
				point_count = i;
				break;
			}
			
			skip(&cur,'\r');
			
			if(skip(&cur,'\n')) {
				*(cur - 1) = 0;
				dim_count = j + 1;
				break;
			}
			
			estimate[i][j][0] = 1.f;
		}
	}
	if(cur != max) {
		fprintf(stderr,"Input file has unexpected format.\n");
		exit(1);
	}
	
	
	u8 new_coord[D] z1;
	
	for(int i = 0; i < point_count; ++i) {
		
		int next = 0;
		int curr;
		int new_node = 0;
		
		for(int j = 0; j < dim_count; ++j) {
			
			new_coord[j] = clamp(0,M,(int) (point[i][j] * delta_m));
			
			curr = next;
			next = index[j][curr][new_coord[j]];
			
			if(!next) {
				
				if(!new_node) {
					new_node = (++node_count) - 1;
				}
				
				next = new_node;
				index[j][curr][new_coord[j]] = next;
			}
		}
		point_node[i] = next;
		memcpy(coord[next],new_coord,sizeof(coord[0]));
		++density_scaled[next];
		
		
		for(int j = 0; j < dim_count; ++j) {
			
			if(estimate[next][j][0] > point[i][j]) {
				estimate[next][j][0] = point[i][j];
			}
			
			if(estimate[next][j][1] < point[i][j]) {
				estimate[next][j][1] = point[i][j];
			}
		}
	}
	
	printf("Stage 1:       #dimensions(d)=%-3d  #points(n)=%-7d  #nodes=%d\n",dim_count,point_count,node_count - 1);
	printf("\n");
}


//  stage 2: process  ----------------------------------------------------------


int adjacent_fast(int a,int b) {
	
	
	// estimate if nodes are possibly adjacent
	//
	//
	// return 1 if maximum_distance <= tau + 0...2*delta,
	// return 0 otherwise
	//
	//
	// equivalent:
	//
	// maximum_distance <= tau + 0...2*delta
	// maximum_node_distance <= tau_factor
	
	
	for(int i = 0; i < dim_count; ++i) {
		
		if(coord[a][i] + tau_factor < coord[b][i]
		|| coord[b][i] + tau_factor < coord[a][i]) {
			return 0;
		}
	}
	
	return 1;
}


int adjacent_exact(int a,int b) {
	
	
	// determine if nodes are adjacent
	//
	//
	// return 1 if maximum_distance <= tau,
	// return 0 otherwise
	
	
	float dist = 0.f,diff;
	
	for(int k = 0; k < dim_count; ++k) {
		
		diff = estimate[a][k][0] - estimate[b][k][1];
		if(dist < diff) {
			dist = diff;
		}
		
		diff = estimate[b][k][0] - estimate[a][k][1];
		if(dist < diff) {
			dist = diff;
		}
	}
	
	if(dist > tau) {
		return 0;
	}
	
	return 1;
}


void* worker_2(void* ptr) {
	
	
	// variants 2,3: cycle nodes
	//
	// iterate nodes, determine if adjacent
	//
	//
	// cost = n^2
	
	
	int t = (int) (long int) ptr;
	for(int i = task[t]; i < task[t + 1]; ++i) {
		
		
		int j = node_count - 1;
		__int128_t *ptr1 = (__int128_t*) coord[i];
		__int128_t *ptr2 = (__int128_t*) coord[1];
		__int128_t xmm0;
		__int128_t xmm1 = *ptr1;
		__int128_t xmm2 = tau_factor;
		__int128_t xmm3 = tau_factor * 2;
		
		
		__asm__ (


"               pxor        %[xmm0],%[xmm0]             ;"
"               pshufb      %[xmm0],%[xmm2]             ;"
"               pshufb      %[xmm0],%[xmm3]             ;"
"               psubb       %[xmm2],%[xmm1]             ;"


			: [xmm0]"=x"(xmm0), [xmm1]"+x"(xmm1), [xmm2]"+x"(xmm2), [xmm3]"+x"(xmm3)
			: 
			: "cc"
		);
		
		
		for(int j2 = 1; j2 < node_count && j > 0; ++j2, --j) {
			
			
			__asm__ (


".align 16                                              ;"
"continue:      movdqa      (%%rdx),%[xmm0]             ;"
"               add         %%rbx,%%rdx                 ;"
"               psubb       %[xmm1],%[xmm0]             ;"
"               psubusb     %[xmm3],%[xmm0]             ;"
"               ptest       %[xmm0],%[xmm0]             ;"
"                                                       ;"
"               jz break                                ;"
"                                                       ;"
"               dec         %%rcx                       ;"
"                                                       ;"
"               jnz continue                            ;"
"                                                       ;"
".align 16                                              ;"
"break:                                                 ;"


				: "+c"(j), "+d"(ptr2), [xmm0]"=&x"(xmm0)
				: "b"(16), [xmm1]"x"(xmm1), [xmm3]"x"(xmm3)
				: "cc"
			);
			
			
			if(j) {
				
				int a = i;
				int b = node_count - j;
				
				if(variant == 2 && !adjacent_exact(a,b)) {
					continue;
				}
				
				if(a != b) {
					edge[a][ edge_count[a]++ ] = b;
					
					if(edge_count[a] == E) {
						fprintf(stderr,"Maximum number of edges per node exceeded.\n");
						exit(1);
					}
				}
			}
		}
	}
	
	return 0;
}


void iterate_recursive(int a,u8* coord_a,u16* curr,int depth) {
	
	
	// iterate adjacent coords, recursive implementation
	
	
	for(int i = -tau_factor; i <= tau_factor; ++i) {
		
		int select = coord_a[depth] + i;
		
		if(select < 0 || select >= M) {
			continue;
		}
		
		
		int next = index[depth][curr[depth]][select];
		
		if(!next) {
			continue;
		}
		
		
		if(depth < dim_count - 1) {
			
			curr[depth + 1] = next;
			iterate_recursive(a,coord_a,curr,depth + 1);
			
		}
		else {
			
			int b = next;
			
			if(!adjacent_exact(a,b)) {
				continue;
			}
			
			if(a != b) {
				edge[a][ edge_count[a]++ ] = b;
				
				if(edge_count[a] == E) {
					fprintf(stderr,"Maximum number of edges per node exceeded.\n");
					exit(1);
				}
			}
		}
	}
}


void* worker_4(void* ptr) {
	
	
	// variant 4: cycle coords
	//
	// iterate adjacent coords, look up nodes
	//
	//
	// cost = n*(t^(d+1)-1)/(t-1) < 2*n*t^d,
	// t = 2*tau_factor+1, skip empty subtrees
	
	
	int t = (int) (long int) ptr;
	for(int i = task[t]; i < task[t + 1]; ++i) {
		
		u8 coord_a[D];
		memcpy(coord_a,coord[i],sizeof(coord[0]));
		
		u16 curr[D];
		curr[0] = 0;
		
		iterate_recursive(i,coord_a,curr,0);
		
	}
	
	return 0;
}


void push(u16 a) {
	stack[ stack_count++ ] = a;
}


u16 pop() {
	return stack[ --stack_count ];
}


void search_inline() {
	
	
	// perform depth-first search to find connected components
	//
	// determine adjacent nodes in one pass
	//
	//
	// for nodes with density > rho, determine adjacent nodes
	// consider as connected if a path of adjacent nodes exists
	//
	//
	// equivalent:
	//
	// density >= rho
	// density*scale >= rho*scale
	//
	// where:
	//
	// scale = n/delta_m^d
	
	
	for(int i = 1; i < node_count; ++i) {
		
		if(!node_component[i]
		&& density_scaled[i] >= rho_scaled) {
			
			node_component[i] = (++component_count) - 1;
			push(i);
			
			for(int j2 = 1; j2 < node_count && stack_count; ++j2) {
				
				int j = pop();
				
				for(int k = i + 1; k < node_count; ++k) {
					
					if(!node_component[k]
					&& density_scaled[k] >= rho_scaled
					&& adjacent_fast(j,k)
					&& adjacent_exact(j,k)) {
						
						node_component[k] = component_count - 1;
						push(k);
						
					}
				}
			}
		}
	}
}


void search_collect() {
	
	
	// perform depth-first search to find connected components
	//
	// use previously generated edge list to determine adjacent nodes
	// see also description of search_inline()
	
	
	for(int i = 1; i < node_count; ++i) {
		
		if(!node_component[i]
		&& density_scaled[i] >= rho_scaled) {
			
			node_component[i] = (++component_count) - 1;
			push(i);
			
			for(int j2 = 1; j2 < node_count && stack_count; ++j2) {
				
				int j = pop();
				
				for(int k2 = 0; k2 < edge_count[j]; ++k2) {
					
					int k = edge[j][k2];
					
					if(!node_component[k]
					&& density_scaled[k] >= rho_scaled) {
						
						node_component[k] = component_count - 1;
						push(k);
						
					}
				}
			}
		}
	}
}


void prefetch(void* ptr,int size) {
	
#ifdef PREFETCH
	
	u8 *cur = (u8*) ptr;
	for(int i = 0; i < size; i += 64,cur += 64) {
		__builtin_prefetch(cur);
	}
	
#endif
	
}


void stage2_process() {
	
	
	// stage 2: process data
	//
	// apply algorithms for density-based cluster analysis
	//
	//
	// variant 0:
	//
	// perform no processing, state no clusters found
	//
	//
	// variant 1:
	//
	// for each rho,
	// determine adjacent nodes and connected components
	//
	//
	// variants 2,3,4:
	//
	// generate list of edges to adjacent nodes
	//
	// for each rho,
	// use edge list to determine connected components
	
	
	if(variant == 0) {  // 0
		return;
	}
	
	if(variant != 1) {  // 2,3,4
		
		if(variant != 4) {  // 2,3
			prefetch(coord,sizeof(coord[0]) * node_count);
		}
		
		task[0] = 1;
		
		for(int t = 0; t < T; ++t) {
			
			task[t + 1] = task[t] + ((node_count - 2) / T) + 1;
			
			if(task[t + 1] > node_count || t == T - 1) {
				task[t + 1] = node_count;
			}
			
			if(variant != 4) {  // 2,3
				pthread_create(&thread[t],NULL,worker_2,(void*) (long int) t);
			}
			else {              // 4
				pthread_create(&thread[t],NULL,worker_4,(void*) (long int) t);
			}
		}
		
		for(int t = 0; t < T; ++t) {
			pthread_join(thread[t],NULL);
		}
	}
	
	
	for(int i = 0; i < point_count && cluster_count == 1; ++i) {
		
		component_count = 1;
		cluster_count = 0;
		memset(node_component,0,sizeof(node_component));
		memset(component_cluster,0,sizeof(component_cluster));
		
		
		// iterate over rho = i/(n*(2*delta)^d) = i/scale
		// scale = n/delta_m^d
		
		rho_scaled = i;
		
		
		if(variant == 1) {
			search_inline();
		}
		else {  // 2,3,4
			search_collect();
		}
		
		
		for(int j = 1; j < node_count; ++j) {
			
			if(node_component[j]) {
				
				
				// equivalent:
				//
				// density >= rho + 2*eps
				// density*scale >= rho*scale + 2*eps*scale
				// density*scale >= rho*scale + eps_factor
				//
				// where:
				//
				// scale = n/delta_m^d
				//
				// max_density = delta_m^d
				// max_density*scale/n = 1
				//
				// 2*eps = sqrt(max_density/(n^2*(2*delta)^d))*eps_factor
				// 2*eps*scale = sqrt(max_density*scale/n)*eps_factor
				// 2*eps*scale = eps_factor
				
				
				if(density_scaled[j] >= rho_scaled + eps_factor
				&& !component_cluster[node_component[j]]) {
					component_cluster[node_component[j]] = ++cluster_count;
				}
			}
		}
	}
	
	printf("Stage 2:       rho*scale=%-8d  #components=%-6d  #clusters=%d\n",rho_scaled,component_count - 1,cluster_count);
	
	if(cluster_count) {
		int node_stats[3] z1;
		for(int i = 0; i < 3; ++i) {
			for(int j = 0; j < point_count; ++j) {
				if(component_cluster[node_component[point_node[j]]] == i) {
					++node_stats[i];
				}
			}
		}
		
		printf("               #cluster1=%-8d  #cluster2=%-8d  #none=%d\n",node_stats[1],node_stats[2],node_stats[0]);
	}
	
	printf("\n");
}


//  stage 3: write  ------------------------------------------------------------


void stage3_write() {
	
	
	// stage 3: write output file
	//
	//
	// prepend cluster id to lines from input
	// if not in cluster, prepend "0"
	// if no clusters found, prepend "1"
	//
	//
	// look up cluster id:
	// point -> node -> component -> cluster
	
	
	FILE *fd2 = fopen(file2,"w");
	
	if(!fd2) {
		fprintf(stderr,"Failed to create specified file.\n");
		exit(1);
	}
	
	char *cur = (char*) buffer;
	char *max = (char*) (buffer + sizeof(buffer));
	
	for(int i = 0; i < point_count && cur < max; ++i) {
		
		int point_cluster = 1;
		
		if(cluster_count) {
			point_cluster = component_cluster[node_component[point_node[i]]];
		}
		
		fprintf(fd2,"%d, ",point_cluster);
		fputs(cur,fd2);
		fputc('\n',fd2);
		
		cur += strlen(cur) + 1;
	}
	
	fclose(fd2);
}


//  program entry point  -------------------------------------------------------


void cluster() {
	
	stage1_read();
	
	stage2_process();
	
	stage3_write();
}


void print_help() {
	fprintf(stderr,"\n");
	fprintf(stderr,"variant=0...4  delta_m=1...64  tau_factor=1,2  eps_factor=0...n  file1,file2=...\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"Read input from file1, write result to file2.\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"delta=1/(2*delta_m)    tau=2*delta*tau_factor\n");
	fprintf(stderr,"2*eps=sqrt(max_density/(n^2*(2*delta)^d))*eps_factor\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"variant 0:   nothing\n");
	fprintf(stderr,"variant 1:   sequential   cycle nodes\n");
	fprintf(stderr,"variant 2:   parallel     cycle nodes    accelerated\n");
	fprintf(stderr,"variant 3:   variant 2 + fast approximation\n");
	fprintf(stderr,"variant 4:   parallel     cycle coords\n");
}


int main(int argc,char* argv[]) {
	
	if(argc != 7) {
		fprintf(stderr,"Exactly 6 command line arguments expected.\n");
		print_help();
		exit(1);
	}
	
	variant = atoi(argv[1]);
	
	delta_m = atoi(argv[2]);
	tau_factor = atoi(argv[3]);
	eps_factor = atoi(argv[4]);
	
	file1 = argv[5];
	file2 = argv[6];
	
	if(variant < 0
	|| variant > 4
	|| delta_m < 1
	|| delta_m > M
	|| tau_factor < 1
	|| tau_factor > 2
	|| eps_factor < 0
	|| eps_factor > N) {
		fprintf(stderr,"Invalid command line argument.\n");
		print_help();
		exit(1);
	}
	
	delta = 0.5f / delta_m;
	tau = 2.f * delta * tau_factor;
	
	printf("variant=%-13d  delta_m=%-6d  tau_factor=%-3d  eps_factor=%d\n",variant,delta_m,tau_factor,eps_factor);
	printf("scale=n/delta_m^d      delta=%-8.4f  tau=%-10.4f  2*eps*scale=%d\n",delta,tau,eps_factor);
	printf("\n");
	printf("file1=%s\n",file1);
	printf("file2=%s\n",file2);
	printf("\n");
	
	cluster();
	
	printf("Success.\n");
	
	exit(0);
}


//  end of code  ---------------------------------------------------------------

