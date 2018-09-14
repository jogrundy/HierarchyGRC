//
//  hierarchy.c
//  expt1
//
//  Created by Joanna Houghton on 07/08/2018.
//  Copyright © 2018 JoSoft. All rights reserved.
//  following  E. Mones, L. Vicsek, and T. Vicsek. Hierarchy measure for complex networks. PLOS ONE, 7(3):1–10, 03 2012.


#include "hierarchy.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utilities.h"
#include <stdbool.h>
#include "matrix_utilities.h"

struct dfs_node{
    int node_n;
    struct dfs_node *prev;
    struct dfs_node *next;
};

void dfs_push(struct dfs_node **ptail, int node_n){
    
    struct dfs_node *new = malloc(sizeof(struct dfs_node));
    new->node_n =node_n;
    if (*ptail) {
        new->prev = *ptail;
        new->next = (*ptail)->next;
        new->next->prev = new;
        new->prev->next = new;
    } else {
        new->next = new;
        new->prev = new;
    }
    *ptail = new;
}

struct dfs_node  * dfs_pop(struct dfs_node **ptail){
    assert(*ptail);
    struct dfs_node *old = *ptail;
    if (old->next == old){
        *ptail = NULL;
    } else {
        (*ptail)->prev->next = (*ptail)->next;
        (*ptail)->next->prev = (*ptail)->prev;
        *ptail = (*ptail)->prev;
    }
    return old;
}


struct bfs_node{
    int node_n;
    int dist;
    double weight;
    struct bfs_node *prev;
    struct bfs_node *next;

};

void enqueue(struct bfs_node **ptail, int node_n, int dist, double weight){
    struct bfs_node *new = malloc(sizeof(struct bfs_node));
    new->node_n = node_n;
    new->dist = dist;
    new->weight = weight;
    if (*ptail) {
        new->next = *ptail;
        new->prev = (*ptail)->prev;
        (*ptail)->prev = new;
        new->prev->next = new;
    } else {
        new->next = new;
        new->prev = new;
        *ptail = new;
    }

}

struct bfs_node  * dequeue(struct bfs_node **ptail){
    assert(*ptail);
    struct bfs_node *old = *ptail;
    if (old->next == old){
        *ptail = NULL;
    } else {
        (*ptail)->prev->next = (*ptail)->next;
        (*ptail)->next->prev = (*ptail)->prev;
        *ptail = (*ptail)->next;
    }
    return old;
}



double get_H(double *matrix, const int N, bool directed, bool weighted){
    // takes in a matrix representing a network
    // assumed to be directed.
    // using the adjacency network removes negative weight problem
    // outputs a score based on GRC, global reaching centrality
    // using local reaching centrality Cr(i) = 1/(N-1) sum (1/d(i,j))
    // for unweighted directed networks Cr(i) = 1/(N-1) * sum (number of reachable nodes)
    printf("network has %i edges and %i nodes\n", get_num_edges(matrix, N*N), N);
    
    double *local_centralities;
    local_centralities = (double *)malloc(N*sizeof(double));
    int i;
    

    if (weighted){
        double *norm_matrix;
        norm_matrix = (double *)malloc(N*N*sizeof(double));
        normalise_weight_matrix(matrix, N*N, norm_matrix);
        for (i=0;i<N;i++){
            local_centralities[i]=local_reaching_centrality_weighted(norm_matrix, N, i, directed);
        }
        free(norm_matrix);

    } else {
        double *adj_mat;
        adj_mat = (double *)malloc(N*N*sizeof(double));
        double *adj_mat_T;
        adj_mat_T = (double *)malloc(N*N*sizeof(double));
        double eps = array_stddev(matrix, N*N)*0.5;
        make_double_adj_mat(matrix, adj_mat, N*N, eps, N, directed);
        matrix_transpose_d(adj_mat, adj_mat_T, N);
        for (i=0;i<N;i++){
            local_centralities[i]=local_reaching_centrality_unweighted(adj_mat, N, i, directed);
        }
        free(adj_mat);
        free(adj_mat_T);
    }
    double sum_diff=0;
    int max_i = argmax(local_centralities, N);
    for (i=0;i<N;i++){
        sum_diff += (local_centralities[max_i] - local_centralities[i]);
    }
    free(local_centralities);
//    double total = total_absolute_weight_array(matrix, N*N);
    return sum_diff/(double)(N-1);
}

double local_reaching_centrality_weighted(double *matrix, const int N, int node_i, bool directed){
//    using local reaching centrality Cr(i) = 1/(N-1) sum (1/d(i,j))
    if (directed){
         // weighted and directed
            double path_weight_sum = sum_over_paths_i(matrix, N, node_i, directed);
            double C_ri = (path_weight_sum)/(double)(N-1);
//        printf("C_ri for node %i is %0.03f\n", node_i, C_ri);
            return C_ri;
    } else { // should work for unweighted and undirected anyway.. should generalise fine.
        // todo: undirected and weighted
    }
    return -1;
}
double local_reaching_centrality_unweighted(double *matrix, const int N, int node_i, bool weighted){
    //    using local reaching centrality Cr(i) = 1/(N-1) sum (1/d(i,j))
    if (weighted){
        // weighted and directed
        double path_weight_sum = sum_over_paths_i(matrix, N, node_i, weighted);
        double C_ri = (path_weight_sum)/(double)(N-1);
        //        printf("C_ri for node %i is %0.03f\n", node_i, C_ri);
        return C_ri;
            //        int sum_over_all_paths = 0;
            //        double C_ri = (1/(double)(N-1))*((double)path_sum);
        } else {
            double path_weight_sum = sum_over_paths_i(matrix, N, node_i, weighted);
            double C_ri = (path_weight_sum)/(double)(N-1);
            //        printf("C_ri for node %i is %0.03f\n", node_i, C_ri);
            return C_ri;

            
        }
    return -1;
    }

double sum_over_paths_i(double *matrix, const int N, int node_i, bool weighted){
    // for loop over all nodes to calculate shortest paths to node_i
    // must first calculate the shortest path between node i and j
    // then calculate it for all j
    int j;
    double weight=0;
    double path_weight_sum = 0;
    if (weighted){
        for (j=0;j<N;j++){
//            printf("node_i = %i, node_j = %i\n", node_i, j);
            int dist_ij = shortest_path_bi(matrix, N, node_i, j, &weight);
            if (dist_ij>0){
//                printf("dist = %i, weight = %0.03f\n", dist_ij, weight);
                path_weight_sum += weight / (double )dist_ij;
//                printf("path_weight_sum += %0.03f so is now: %0.03f\n", weight / (double )dist_ij, path_weight_sum);
            }
        }
    } else {
        for (j=0;j<N;j++){
            //            printf("node_i = %i, node_j = %i\n", node_i, j);
            int dist_ij = shortest_path_bi(matrix, N, node_i, j, &weight);
            if (dist_ij>0){
                //                printf("dist = %i, weight = %0.03f\n", dist_ij, weight);
                path_weight_sum += 1;
            }
        }
    }
//    printf("path_weight_sum for node %i = %0.03f\n",  node_i, path_weight_sum);
//    double total = total_absolute_weight_array(matrix, N*N);
    return path_weight_sum;
}

double proportion_reachable(int *matrix, const int N, int node_i){
    // returns proportion reachable from a particular node node_i
    int j, total=0;
    for (j=0;j<N;j++){
        if (reachable(matrix, N, node_i, j)==1){
            total ++;
        }
    }
    return (double)total/(double )(N-1);
}

//int shortest_path(double *matrix, const int N, int node_i, int node_j, int *path, bool directed){
//    // does bfs over nodes to find shortest path from i to j
//    // uses bfs node to create data structure - last in last out
//    // will only work for small simple networks, memory will become problem quickly.
//    // would it be better just to say that the node has been found?
//    // as I'm not using the path distance at the moment..
//
//    if (node_i == node_j) {
//        path[0] = node_i;
//        return 0; // if start and end nodes are the same, then path length is 0
//    }
//    struct bfs_node *headtail = NULL; // start queue
//    int *visited;
//    int dist = 0;
//    visited = (int *)malloc(N*sizeof(int));
//    int n_visited = 1;
//    enqueue(&headtail, node_i, dist, 0);
//
//    int *out_nodes;
//    out_nodes = (int *)malloc(N*sizeof(int));
//
//    int count = 1;
//    while (headtail) {
//        struct bfs_node *new = dequeue(&headtail);
//
//        int node_n = new->node_n;
//        int dist = new->dist;
//        free(new);
//        count --;
//        int i;
//        // find where you can go from this node
//
//        int n_out = find_out_nodes_unweighted(matrix, N, node_n, out_nodes);
//        //explore each possible path
//
//        for (i=0;i<n_out;i++){
//            int next_node_n = out_nodes[i];
//            if (next_node_n == node_j){ // found it!
//                free(out_nodes);
//                free(visited);
//                clear_queue(&headtail);
////                printf("count = %i\n", count);
//                return dist;
//            } else if (is_in(visited, next_node_n, dist+1)){
//                // been there before
//                continue;
//            } else {
//                enqueue(&headtail, node_n, dist+1, 0);
//                count ++;
//                visited[n_visited] = node_n;
//                n_visited ++;
//            }
//
//        }
//
//    }
//    free(out_nodes);
//    free(visited);
//    // no path is possible, so path length infinite?? better to only use in connected network bits??
//    return -1;
//}

int shortest_path_bi(double *matrix, const int N, int node_i, int node_j, double *weight){
    // does bfs over nodes to find shortest path from i to j
    // uses bfs node to create data structure - last in last out
    // modified to use bidirectional search
    
//    printf("node_i=%i, node_j=%i\n", node_i, node_j);
    if (node_i == node_j) {
        weight = 0;
//        printf("found node! dist=%i\n\n", 0);
        return 0; // if start and end nodes are the same, then path length is 0
    }
    struct bfs_node *start_Q = NULL; // start queue
    struct bfs_node *goal_Q = NULL; // start queue

    enqueue(&start_Q, node_i, 0, 0.00); // put first node in queue
    enqueue(&goal_Q, node_j, 0, 0.00); // put first node in queue
    
    int start_count = 1; // counting up num nodes in the queue
    int goal_count = 1;
    
     // list of nodes visited from start
    int *start_list;
    start_list = (int *)malloc(N*sizeof(int));
    start_list[0] = node_i;
    int n_start_list = 1; // counting number in visited list
    
    // list of nodes visited from goal
    int *goal_list;
    goal_list = (int *)malloc(N*sizeof(int));
    goal_list[0]= node_j; // list of nodes visited from end
    int n_goal_list = 1;
    
    //allocate memory for in nodes and out nodes lists
    int *in_nodes;
    in_nodes = (int *)malloc(N*sizeof(int));
    int *out_nodes;
    out_nodes = (int *)malloc(N*sizeof(int));
    

    while (start_Q || goal_Q) {
        // look at shortest queue first
        // need visited list for both queues so i can check if theyve met in the middle
        struct bfs_node *new;
        int node_n;
        int dist;
        double node_weight;

        if (start_Q && (n_start_list<=n_goal_list || !goal_Q)){ // get new node to expand,
//            printf("start_Q first\n");
            new = dequeue(&start_Q);
            node_n = new->node_n;
            dist = new->dist;
            node_weight = new->weight;
            free(new);
            start_count --;
//            printf("expanding node %i\n", node_n);
            // find where you can go from this node
            int n_out = find_out_nodes_weighted(matrix, N, node_n, out_nodes);
            int i;
            for (i=0;i<n_out;i++){
                int next_node_n = out_nodes[i];
//                printf("can get to node %i from node %i\n", next_node_n, node_n);
                
                // goal checking..
                if (is_in(goal_list, next_node_n, n_goal_list)){ // found a path!
                    bool found = false;
                    
                    int g_node_n=0;
                    int g_dist=0;
                    double g_weight=0;
                    while(!found){
                        struct bfs_node *next_one = dequeue(&goal_Q);
                        if (next_one->node_n == node_n ||next_one->node_n == next_node_n || next_one->node_n == node_j) {// looking for the member of the goal queue that matches the node in the start queue
                            g_node_n = next_one->node_n;
                            g_dist = next_one->dist;
                            g_weight = next_one->weight;
                            free(next_one);
                            found = true;
                        } else {
                            free(next_one);
                        }

                    }
                    dist += g_dist;
                    dist +=1;

                    node_weight += g_weight;
                    node_weight += give_weight_i_j(matrix, N, node_n, g_node_n);

                    *weight = node_weight;
//                    printf("found node! dist=%i path weight = %0.03f\n\n", dist, *weight);
//                    if (dist != *weight){
//                        printf("oops something has gone wrong here\n");
//                    }
                    
                    clear_queue(&start_Q);
                    clear_queue(&goal_Q);
                    free(start_list);
                    free(goal_list);
                    free(in_nodes);
                    free(out_nodes);
                    return dist; // 1 for found
                    
                } else if (is_in(start_list, next_node_n, n_start_list)){
                    // been there before
                    continue;
                    
                } else {
                    // add to queue, add to visited start_list, increase start_count
                    double weight_ij = give_weight_i_j(matrix, N, node_n, next_node_n);
                    enqueue(&start_Q, next_node_n, dist+1, node_weight+weight_ij);
//                    printf("enqueueing node %i, with dist=%i and total path weight so far %0.03f\n", next_node_n, dist+1, node_weight+weight_ij);
                    start_count++;
                    start_list[n_start_list] = next_node_n;
                    n_start_list++;
                }
            }
        } else {
            // goal is shortest queue get new node out of queue, get new_node number, free new
//            printf("goal Q first\n");
            new = dequeue(&goal_Q);
            node_n = new->node_n;
            dist = new->dist;
            node_weight = new->weight;
            free(new);
//            printf("expanding node %i\n", node_n);

            goal_count --;
            
            // find where you can go from this node

            int n_in = find_in_nodes_weighted(matrix, N, node_n, in_nodes);
            int i;
            for (i=0;i<n_in;i++){
                //explore each possible path
                
                int next_node_n = in_nodes[i];
//                printf("can get to node %i from node %i\n", node_n, next_node_n);
                
                if (is_in(start_list, next_node_n, n_start_list)){ // found a path!
                    bool found = false;
                    int s_node_n=0;
                    int s_dist=0;
                    double s_weight=0;
                    while(!found){
                        struct bfs_node *next_one = dequeue(&start_Q);
                        if (next_one->node_n == node_n ||next_one->node_n == next_node_n || next_one->node_n == node_i)  {// looking for the member of the start queue that matches the node in the start queue
                            s_node_n = next_one->node_n;
                            s_dist = next_one->dist;
                            s_weight = next_one->weight;
                            free(next_one);
                            found = true;
                        } else {
                            free(next_one);
                        }
                    }
                    dist += s_dist;
                    dist +=1;
                    
                    node_weight += s_weight;
                    node_weight += give_weight_i_j(matrix, N, s_node_n, node_n);
                    *weight = node_weight;
//                    printf("found node! dist=%i path weight = %0.03f\n\n", dist, node_weight);
//                    if (dist != *weight){
//                        printf("oops something has gone wrong here\n");
//                    }
                    clear_queue(&start_Q);
                    clear_queue(&goal_Q);
                    free(start_list);
                    free(goal_list);
                    free(in_nodes);
                    free(out_nodes);
                    return dist; // 1 for found

                } else if (is_in(goal_list, next_node_n, n_goal_list)){
                    // been there before
                    continue;
                } else {
                    double weight_ij = give_weight_i_j(matrix, N, next_node_n, node_n);
                    enqueue(&goal_Q, next_node_n, dist+1, node_weight+weight_ij);
                    goal_count ++;
//                    printf("enqueueing node %i, with dist=%i and total path weight so far %0.03f\n", next_node_n, dist+1, node_weight+weight_ij);
                    // update visited list
                    goal_list[n_goal_list] = next_node_n;
                    n_goal_list ++;
                }
            }
        }
    }
    // no path is possible
    free(start_list);
    free(goal_list);
    free(in_nodes);
    free(out_nodes);
//    printf("node is unreachable\n\n");
    return 0;
}

void clear_queue(struct bfs_node **headtail){
    // to make sure all memory is freed when the target is found
    while(*headtail){
        struct bfs_node *surplus = dequeue(headtail);
        free(surplus);
    }
}

void empty_heap(struct dfs_node **headtail){
    // to make sure all memory is freed when the target is found
    while(*headtail){
        struct dfs_node *surplus = dfs_pop(headtail);
        free(surplus);
    }
}



int find_in_nodes_weighted(double *matrix, const int N, int node_i, int *in_nodes){
//    takes in matrix and node number, returns number of in connections
//    in_nodes carries the node numbers of these nodes
    // uses weighted connection matrix.
    int j, n_in=0;
    memset(in_nodes, 0, N*sizeof(int)); // could be N big, but don't know yet.
    for (j=0;j<N;j++){
        int ind = node_i*N + j;
        if (matrix[ind]>0){
            in_nodes[n_in] = j;
            n_in ++;
        }
        
        } // basically just returns a node list and a number nodes n_in = N.
    // depends on the sparsity of the matrix.
    return n_in;
}

int find_in_nodes_unweighted(int *matrix, const int N, int node_i, int *in_nodes){
    //    takes in matrix and node number, returns number of in connections
    //    in_nodes carries the node numbers of these nodes
    // uses adjacency matrix
    int j, n_in=0;
    memset(in_nodes, 0, N*sizeof(int)); // could be N big, but don't know yet.
    for (j=0;j<N;j++){
        int ind = node_i*N + j;
        if (matrix[ind] > 0) {
            in_nodes[n_in] = j;
            n_in ++;
        }
    }
    return n_in;
}
int find_out_nodes_unweighted(int *matrix, const int N, int node_i, int *out_nodes){
    //    takes in matrix and node number, returns number of out connections
    //    out_nodes carries the node numbers of these nodes
    // using adjacency matrix
    int i, n_out=0;
    memset(out_nodes, 0, N*sizeof(int)); // could be N big, but don't know yet.
    for (i=0;i<N;i++){
        int ind = i*N + node_i;
//        printf("matrix at ind %i = %f\n", ind, matrix[ind]);
        if (matrix[ind] > 0) {
            out_nodes[n_out] = i;
            n_out ++;
        }
    }
    return n_out;
}

int find_out_nodes_weighted(double *matrix, const int N, int node_i, int *out_nodes){
    //    takes in matrix and node number, returns number of in connections
    //    in_nodes carries the node numbers of these nodes
    // uses adjacency matrix
    int j, n_out=0;
    memset(out_nodes, 0, N*sizeof(int)); // could be N big, but don't know yet.
    for (j=0;j<N;j++){
        int ind = j*N+node_i;
        if (matrix[ind]>0){
            out_nodes[n_out] = j;
            n_out ++;
        }
    }
    return n_out;
}

int reachable (int *matrix, int N, int node_i, int node_j){
    // does dfs over nodes to find shortest path from i to j
    // uses dfs node to create data structure - last in first out
    // hopefully better for memory than BFS search, and i only need to know if it is reachable

//    printf("node_i = %i, node_j = %i\n", node_i, node_j);
    if (node_i == node_j) {
        return 1; // return 1 for reachable and -1 for not reachable
    }
    struct dfs_node *headtail = NULL; // start queue
    dfs_push(&headtail, node_i); // put first node in queue

    int count = 1;
//    printf("starting at node %i, going to node %i\n", node_i, node_j);
    // visited list
    int *visited;
    int n_visited=1;
    visited = (int *)malloc(N*sizeof(int)); // cannot be bigger than N
    visited[0] = node_i;
    
    int *out_nodes;
    out_nodes = (int *)malloc(N*sizeof(int)); // cannot be bigger than N
    
    while (headtail) {
        
        struct dfs_node *new = dfs_pop(&headtail);
        int node_n = new->node_n;
        free(new); // everything that is popped is freed
        count --;

        // find where you can go from this node
        int n_out = find_out_nodes_unweighted(matrix, N, node_n, out_nodes);
        
        // explore each possible path
        int i;
        for (i=0;i<n_out;i++){
            int next_node_n = out_nodes[i];
            if (next_node_n == node_j){ // found it!
                //free all allocated memory
                free(out_nodes);
                free(visited);
                empty_heap(&headtail);
                return 1;
                
            } else if (is_in(visited, next_node_n, n_visited)){
                // been there before
                continue;
            } else {
                
                dfs_push(&headtail, next_node_n);
                count ++;
                // add node to visited list, increase visited count
                visited[n_visited] = node_n;
                n_visited ++;
//                printf("%i\n", node_n);
            }
        }
    }
    // no path is possible, so path length infinite?? better to only use in connected network bits??
    free(out_nodes);
    free(visited);
    return 0;
}

void add_element_to_int_array(int *array, int element, int len){
    // len is length of array before element is added
    array[len] = element;
}

double give_weight_i_j(double *matrix, int N, int node_i, int node_j){
    int ind = node_j*N + node_i;
    return (double )matrix[ind];
}

void normalise_weight_matrix(double *matrix, int len, double *norm_matrix){
    
    // interpret as distance then normalise
    double *inv_matrix;
    inv_matrix = (double *)malloc(len*sizeof(double));
    invert_all_elements(matrix, len, inv_matrix); // now proportional to distance
    
    double total_dist = total_absolute_weight_array(inv_matrix, len); //
    int num_edges = get_num_edges(matrix, len);
    double multiple = (double )num_edges / total_dist;
    multiply_all_elements_d(inv_matrix, len, norm_matrix, multiple); // puts result in norm_matrix
    free(inv_matrix);
}

int get_num_edges(double *matrix, int len){
    int i, n_edges=0;
    for (i=0;i<len;i++){
        if (fabs(matrix[i])>0){
            n_edges ++;
        }
    }
    return n_edges;
}



void test_hierarchy(void){
    const int N1 = 10;
    double test1[N1*N1] ={0, 1, 1, 0, 0, 0, 0, 0, 0, 1,
        1, 0, 1, 0, 0, 0, 0, 0, 0, 0,
        1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 1, 0, 0, 0, 1,
        0, 0, 0, 1, 0, 1, 0, 0, 0, 0,
        0, 0, 0, 1, 1, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 1, 1, 1,
        0, 0, 0, 0, 0, 0, 1, 0, 1, 0,
        0, 0, 0, 0, 0, 0, 1, 1, 0, 0,
        1, 0, 0, 1, 0, 0, 1, 0, 0, 0};
    double test2[N1*N1] ={0, 1, 1, 0, 0, 0, 0, 0, 0, 0,
        1, 0, 1, 0, 0, 0, 0, 0, 0, 0,
        1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 1, 0, 0, 0, 0,
        0, 0, 0, 1, 0, 1, 0, 0, 0, 0,
        0, 0, 0, 1, 1, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 1, 1, 0,
        0, 0, 0, 0, 0, 0, 1, 0, 1, 0,
        0, 0, 0, 0, 0, 0, 1, 1, 0, 0,
        1, 0, 0, 1, 0, 0, 1, 0, 0, 0};
    double test3[N1*N1] ={0, 1, 1, 0, 0, 0, 0, 0, 0, 1,
        1, 0, 1, 0, 0, 0, 0, 0, 0, 0,
        1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 1, 0, 0, 0, 1,
        0, 0, 0, 1, 0, 1, 0, 0, 0, 0,
        0, 0, 0, 1, 1, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 1, 1, 1,
        0, 0, 0, 0, 0, 0, 1, 0, 1, 0,
        0, 0, 0, 0, 0, 0, 1, 1, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double test4[N1*N1] ={0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double test5[N1*N1] ={0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 2,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 2,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 2,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 3,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 2,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 2,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 2,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 4,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
    int N_ce = 297;
    double *celegansnn;
    celegansnn = (double *)malloc(N_ce*N_ce*sizeof(double));
    double *celegansnn_T;
    celegansnn_T = (double *)malloc(N_ce*N_ce*sizeof(double));
    char *filename_ce = "/Users/jojo/MSc Stuff/Project/code/celegansneural.txt";
    load_matrix(celegansnn, N_ce*N_ce, filename_ce);
    matrix_transpose_d(celegansnn, celegansnn_T, N_ce);
    
    int N_little = 186;
    double *littlerock;
    double *littlerock_T;
    littlerock = (double *)malloc(N_little*N_little*sizeof(double));
    littlerock_T = (double *)malloc(N_little*N_little*sizeof(double));
    char *filename_little = "/Users/jojo/MSc Stuff/Project/code/littlerock_array.txt";
    load_matrix(littlerock, N_little*N_little, filename_little);
    matrix_transpose_d(littlerock, littlerock_T, N_little);
    
    bool directed = true;
    bool weighted = false;
    double H1 = get_H(test1, N1, directed, weighted);
    printf("test1 H = %0.03f, should be 0.000\n", H1);
    double H2 = get_H(test2, N1, directed, weighted);
    printf("test2 H = %0.03f, should be 0.037\n", H2);
    double H3 = get_H(test3, N1, directed, weighted);
    printf("test3 H = %0.03f, should be 0.778\n", H3);
    double H4 = get_H(test4, N1, directed, weighted);
    printf("test4 H = %0.03f, chain should be 0.556\n", H4);
    double H5 = get_H(test5, N1, directed, weighted);
    printf("test5 H = %0.03f, star should be 1.0\n", H5);
    double H6 = get_H(celegansnn_T, N_ce, directed, weighted);
    printf("test6 H = %0.03f, celegans nn should be 0.133\n", H6);
    double H7 = get_H(littlerock_T, N_little, directed, weighted);
    printf("test7 H = %0.03f, littlerock should be 0.811\n", H7);
    free(celegansnn);
}
