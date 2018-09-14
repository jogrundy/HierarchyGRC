//
//  hierarchy.h
//  expt1
//
//  Created by Joanna Houghton on 07/08/2018.
//  Copyright © 2018 JoSoft. All rights reserved.
//

#ifndef hierarchy_h
#define hierarchy_h

#include <stdio.h>
#include <math.h>
#include <stdbool.h>
struct bfs_node;
struct dfs_node;
void enqueue(struct bfs_node **ptail, int node_n, int dist, double weight);
double get_H(double *matrix, const int n_nodes, bool directed, bool weighted);
double local_reaching_centrality_weighted(double *matrix, const int N, int node_i, bool directed);
double local_reaching_centrality_unweighted(double *matrix, const int N, int node_i, bool directed);
double sum_over_paths_i(double *matrix, const int N, int node_i, bool directed);
int shortest_path(double *matrix, const int N, int node_i, int node_j, int *path, bool directed);
int shortest_path_bi(double *matrix, const int N, int node_i, int node_j, double *weight);
int find_in_nodes_weighted(double *matrix, const int N, int node_i, int *in_nodes);
int find_in_nodes_unweighted(int *matrix, const int N, int node_i, int *in_nodes);
int find_out_nodes_weighted(double *matrix, const int N, int node_i, int *out_nodes);
int find_out_nodes_unweighted(int *matrix, const int N, int node_i, int *out_nodes);

void make_node_queue_add(struct bfs_node **headtail, const int N, int node_i,  int dist, int *visited);
void test_hierarchy(void);
void clear_queue(struct bfs_node **headtail);
void empty_heap(struct dfs_node **headtail);
int reachable (int *matrix, int N, int node_i, int node_j);
double proportion_reachable(int *matrix, const int N, int node_i);
void dfs_push(struct dfs_node **ptail, int node_n);
struct dfs_node  * dfs_pop(struct dfs_node **ptail);
double give_weight_i_j(double *matrix, int N, int node_i, int node_j);
void normalise_weight_matrix(double *matrix, int len, double *norm_matrix);
int get_num_edges(double *matrix, int len);

#endif /* hierarchy_h */
