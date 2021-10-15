#include <iostream>
#include <fstream>
#include<stdio.h>
#include<stdlib.h>
#include<vector>
#include<limits.h>
#include<float.h>
#include <assert.h>
#include<memory.h>
#include<math.h>
#include<set>
#include<algorithm>
#include<utility>
#include <cstdlib>
#include <cmath>
#include <string.h>
#include <time.h>
#include <unistd.h>

#define MAXA 12

int greedy_limit = 100000;
int cca_limit = 50000;
int smooth = 0;

using namespace std;

double time_limit;

class TimeTest
{
	clock_t start, finish;
public:
	void Init()
	{
		start = clock();
	}
	double GetTime()
	{
		finish = clock();
		return 1.0 * (finish - start) / CLOCKS_PER_SEC;
	}
	bool NeedToStop()
	{
		double nowTime = GetTime();
		return nowTime > time_limit;
	}
} time_test;

typedef struct Column
{
	int score;
	bool config;
	bool is_in_c;
} Column;

typedef struct Row
{
	int weight;
} Row;

Column **cs;//set
Row **rs;//element
int m, n;
int element_size, set_size;
int seed;
int density;
int **uncover;
int **index_uncover;
int **subscore;
int **covered;
int *uncov_num;
int *Remember;
int *remember;
int *rtemp;
int *rem_temp;
int *rem_half;
int *col_temp;
int *nncol;
double *tau;
double *eta;
double *te_beta;
int step;
int best_value;
int **ori_row, **ori_col;
int *ori_ncol, *ori_nrow;
int ***row, ***col;
int **ncol, **nrow;
int *best_sol, *k;
double start_time;
int max_k;
int **set_age;
int tabu = -1;
int uncover_num;
double real_time;
int happy;
int sad;
int best_temp;
int greedytag = 0;
double omiga = 0.95;
int randstep = 0;
int cntants = 0;
int tpp = 0;
int uncover_temp = 0;
int tag_init = 1;
int randant = 0;
int valuek = 0;
double finish_time = 0.00;
int randnum = 0;
int beta = 1;
double epslion = 0.005;
int ants = 10;
double delta_tau = 0.;
double taumax, taumin;
int k_half = 0;		//k/2
int remove_score = 0;



void init();
void find_worst_in_c(int);
void remove(int, int);
void find_random_in_vout(int);
void add(int, int);
void reload_weight_score(int);
void reload_age(int);
void CCA_local_search(int);
void CCA_local_searchw(int);
void init_aco();
void init_eta_tebeta();
void add_aoc(int, int);
void reload_alltau();
void reload_alltaupro(int);
void reload_bt(int);
void reload_bs(int);
void maxmin();
void ANT_1();
void ANT_2();
void ANT_3();
void ANT_4();
void ANT_5();
void ANT_6();
void ANT_7();
void ANT_8();
void ANT_greedy();
void ANT_greedy();
void ant_empty_greedy(int t);
void remove_greedy(int t, int c);
void add_aoc_greedy(int t, int c);
void add_greedy(int t,int c);
void ant_search_greedy(int t);
void greedy_ls(int t);
void greedy_worst(int t);
void greedy_in_vout(int t);

void ant_empty(int);
void ant_build(int, int);
void ant_search(int);
void ant_searchw(int);
void init_search(int);
void init_searchw(int);
void MMAS_pro();
void NuAnts();
void update_best_sol();
void free_all();
void free_ori();
void NuAnts_density();
void MMAS_pro_density();




void free_ori()
{
	free(ori_row);
	free(ori_col);
	free(ori_ncol);
	free(ori_nrow);
}

void free_all()
{
	free(row);
	free(col);
	free(ncol);
	free(nrow);
	free(best_sol);
	free(set_age);
	for (int i = 0; i < MAXA; i++)
	{
		free(cs[i]);
		free(rs[i]);

		//        free(uncover[i]);
		//        free(index_uncover[i]);
		delete[] uncover[i];
		delete[] index_uncover[i];

		free(subscore[i]);
		//free(covered[i]);
		//free(Proba[i]);
	}
	//free(best_array);

	free(uncov_num);
	free(Remember);
	free(remember);
	free(rtemp);
	free(rem_temp);
	free(rem_half);
	free(col_temp);
	free(nncol);
	free(tau);
	free(eta);
	free(te_beta);
}


void build_instance(char *file)
{
	int i, j, h, t;
	int joke = 0;
	freopen(file, "r", stdin);
	scanf("%d%d", &m, &n);

	element_size = m + 1;
	set_size = n + 1;

	subscore = (int **)malloc(MAXA * sizeof(int *));
	covered = (int **)malloc(MAXA * sizeof(int *));
	for (i = 0; i < MAXA; i++) subscore[i] = (int *)malloc(set_size * sizeof(int));

	for (i = 0; i < MAXA; i++) covered[i] = (int *)malloc(element_size * sizeof(int));

	for (i = 0; i < MAXA; i++)
	{
		for (j = 0; j < set_size; j++)
		{
			subscore[i][j] = 0;
		}
	}

	for (i = 0; i < MAXA; i++)
	{
		for (j = 0; j < element_size; j++)
		{
			covered[i][j] = 0;
		}
	}

	uncov_num = (int *)malloc(MAXA * sizeof(int));
	for (i = 0; i < MAXA; i++) uncov_num[i] = 0;

	Remember = (int *)malloc(set_size * sizeof(int));
	remember = (int *)malloc(set_size * sizeof(int));
	rtemp = (int *)malloc(set_size * sizeof(int));
	rem_temp = (int *)malloc(element_size * sizeof(int));
	rem_half = (int *)malloc(element_size * sizeof(int));
	col_temp = (int *)malloc(set_size * sizeof(int));
	nncol = (int *)malloc(element_size * sizeof(int));
	for (i = 0; i < element_size; i++)
	{
		nncol[i] = 0;
	}
	for (i = 0; i < set_size; i++)
	{
		col_temp[i] = 0;
	}

	tau = (double *)malloc(set_size * sizeof(double));
	eta = (double *)malloc(set_size * sizeof(double));
	te_beta = (double *)malloc(set_size * sizeof(double));

	cs = (Column **)malloc(MAXA * sizeof(Column *));
	rs = (Row **)malloc(MAXA * sizeof(Row *));
	for (i = 0; i < MAXA; i++)
	{
		cs[i] = (Column *)malloc(set_size * sizeof(Column));
		rs[i] = (Row *)malloc(element_size * sizeof(Row));
	}

	uncover = new int*[MAXA];
	index_uncover = new int*[MAXA];
	for (i = 0; i < MAXA; i++)
	{
		uncover[i] = new int[set_size];
		index_uncover[i] = new int[element_size];
	}


	ori_row = (int **)malloc(m * sizeof(int *));
	ori_ncol = (int *)malloc(m * sizeof(int));
	ori_col = (int **)malloc(n * sizeof(int *));
	ori_nrow = (int *)malloc(n * sizeof(int));
	row = (int ***)malloc(MAXA * sizeof(int **));
	ncol = (int **)malloc(MAXA * sizeof(int *));
	col = (int ***)malloc(MAXA * sizeof(int **));
	nrow = (int **)malloc(MAXA * sizeof(int *));
	set_age = (int **)malloc(MAXA * sizeof(int *));
	k = (int *)malloc(n * sizeof(int));
	best_sol = (int *)malloc(n * sizeof(int));
	for (t = 0; t < MAXA; t++)
	{
		set_age[t] = (int *)malloc(n * sizeof(int));
		for (j = 0; j < n; j++)
			set_age[t][j] = 1;
	}
	for (j = 0; j < n; j++)
	{
		scanf("%d", &joke);
	}

	for (t = 0; t < MAXA; t++) {
		for (j = 0; j < n; j++)
		{
			cs[t][j].config = 1;
			cs[t][j].is_in_c = 0;
			cs[t][j].score = 0;
			best_sol[j] = 0;
		}
	}
	for (i = 0; i < m; i++) 
	{
		scanf("%d", &ori_ncol[i]);
		ori_row[i] = (int *)malloc(ori_ncol[i] * sizeof(int));
		for (h = 0; h < ori_ncol[i]; h++)
		{
			scanf("%d", &ori_row[i][h]);
			ori_row[i][h]--;
		}
	}
	for (j = 0; j < n; j++)
		ori_nrow[j] = 0;
	for (i = 0; i < m; i++)
	{
		for (h = 0; h < ori_ncol[i]; h++)
			ori_nrow[ori_row[i][h]]++;
	}
	for (j = 0; j < n; j++)
	{
		ori_col[j] = (int *)malloc(ori_nrow[j] * sizeof(int));
		k[j] = 0;
	}
	for (i = 0; i < m; i++)
	{
		for (h = 0; h < ori_ncol[i]; h++)
		{
			ori_col[ori_row[i][h]][k[ori_row[i][h]]] = i;
			k[ori_row[i][h]]++;
		}
	}


	for (t = 0; t < MAXA; t++)
	{

		ncol[t] = (int *)malloc(m * sizeof(int));
		nrow[t] = (int *)malloc(n * sizeof(int));
		col[t] = (int **)malloc(n * sizeof(int *));
		row[t] = (int **)malloc(m * sizeof(int *));

		for (i = 0; i < m; i++) {
			rs[t][i].weight = 1;
			ncol[t][i] = ori_ncol[i];
			row[t][i] = (int *)malloc(ori_ncol[i] * sizeof(int));
			for (h = 0; h < ori_ncol[i]; h++)
			{
				row[t][i][h] = ori_row[i][h];
			}
		}

		for (j = 0; j < n; j++)
		{
			cs[t][j].score = ori_nrow[j];
			nrow[t][j] = ori_nrow[j];

			col[t][j] = (int *)malloc(ori_nrow[j] * sizeof(int));
			for (h = 0; h < ori_nrow[j]; h++)
			{
				col[t][j][h] = ori_col[j][h];
			}
		}
	}
	int x = m * n;
	int y = 0;
	for (j = 0; j < n; j++)
	{
		y += ori_nrow[j];
	}
	float z = 1.0 *  y / x;
	printf("m: %d; n: %d; m*n: %d; density: %f \n",m,n,x,z);
	if(z < 0.01)
	{
		density = 1;
		cca_limit = 1000;
		printf("density < 0.01! \n");
	}
	else
		density = 0;
	free(k);
	free_ori();
}



void init()
{
	int i, j;
	uncover_num = m;
	for (i = 0; i < ants; i++) {
		uncov_num[i] = m;
		for (j = 0; j < m; j++)
			index_uncover[i][j] = 1;
	}
}

void Four_stage_solution_structure(int t)
{
	int frcount = 0;
	int i, j;
	int max = -1;
	for (j = 0; j < n; j++)
	{
		if (cs[t][j].is_in_c == 1)
			continue;
		if (cs[t][j].score > max)
		{
			max = cs[t][j].score;
			//cmax = j;
			for (i = 0; i < frcount; i++)
			{
				col_temp[i] = 0;
			}
			frcount = 0;
		}
		if (cs[t][j].score == max)
		{
			col_temp[frcount] = j;
			frcount = frcount + 1;
		}
	}
	if (frcount == 1)
	{
		add_aoc(t, col_temp[0]);
		//SD ++;
	}
	if (frcount > 1)
	{
		int the_set;
		int crcount = 0;
		for (i = 0; i < frcount; i++)
		{
			int max_cnt = 0;
			for (j = 0; j < nrow[t][col_temp[i]]; j++)
			{
				if (index_uncover[t][col[t][col_temp[i]][j]] == -1)
					continue;
				if (ncol[t][col[t][col_temp[i]][j]] == 1)
					max_cnt++;
			}
			if (max_cnt > crcount)
			{
				crcount = max_cnt;
				the_set = col_temp[i];
			}
		}
		if (crcount > 0)
		{
			add_aoc(t, the_set);
		}
		for (i = 0; i < m; i++) {
			nncol[i] = 0;
		}
		if (crcount == 0) {

			int trcount = 0;
			int pse_set;
			for (i = 0; i < frcount; i++) {
				for (j = 0; j < nrow[t][col_temp[i]]; j++) {
					if (index_uncover[t][col[t][col_temp[i]][j]] == -1)
						continue;
					nncol[col[t][col_temp[i]][j]] ++;
				}
			}
			for (i = 0; i < frcount; i++) {
				int pse_cnt = 0;
				for (j = 0; j < nrow[t][col_temp[i]]; j++) {
					if (index_uncover[t][col[t][col_temp[i]][j]] == -1)
						continue;
					if (nncol[col[t][col_temp[i]][j]] == 1)//nncol[]
						pse_cnt++;
				}
				if (pse_cnt > trcount)
				{
					trcount = pse_cnt;
					pse_set = col_temp[i];
				}
			}
			if (trcount > 0) {
				add_aoc(t, pse_set);
			}
			if (trcount == 0) {
				int max_set = 0;
				int max_all = 0;
				for (i = 0; i < frcount; i++) {
					if (max_set < nrow[t][col_temp[i]]) {
						max_set = nrow[t][col_temp[i]];
						max_all = col_temp[i];
					}
				}
				add_aoc(t, max_all);
			}
		}
	}
}

void FourlevelGreedy(int t)
{
	int add_count = 0;
	while (uncover_num > 0)
	{
		if (add_count == max_k)
			break;
		Four_stage_solution_structure(t);
		add_count++;
	}
	update_best_sol();
	init_search(t);
}

void FourlevelGreedyw(int t)
{
	int add_count = 0;
	while (uncover_num > 0)
	{
		if (add_count == max_k)
			break;
		Four_stage_solution_structure(t);
		add_count++;
	}
	update_best_sol();
	init_searchw(t);
}


void find_worst_in_c(int t)
{
	int maxc = 0, maxx = INT_MIN;
	int k_temp = max_k;
	int tempk = 0;
	for (int j = 0; j < k_half; j++)
	{
		randnum = rand() % k_temp;
		rem_half[j] = rem_temp[randnum];
		k_temp--;
		tempk = rem_temp[randnum];
		rem_temp[randnum] = rem_temp[k_temp];
		rem_temp[k_temp] = tempk;
	}
	randstep = rand() % k_half;
	for (int j = randstep; j < k_half; j++)
	{
		if (rem_half[j] == tabu)
			continue;
		if (maxx < cs[t][rem_half[j]].score)
		{
			maxc = rem_half[j];
			maxx = cs[t][rem_half[j]].score;
		}
		else if (maxx == cs[t][rem_half[j]].score)
		{
			if (subscore[t][maxc] > subscore[t][rem_half[j]])
				maxc = rem_half[j];
			else if (subscore[t][maxc] == subscore[t][rem_half[j]])
			{
				if (set_age[t][maxc] < set_age[t][rem_half[j]])
					maxc = rem_half[j];
			}
		}
	}
	for (int j = 0; j < randstep; j++)
	{

		if (rem_half[j] == tabu)
			continue;
		if (maxx < cs[t][rem_half[j]].score)
		{
			maxc = rem_half[j];
			maxx = cs[t][rem_half[j]].score;
		}
		else if (maxx == cs[t][rem_half[j]].score)
		{
			if (subscore[t][maxc] > subscore[t][rem_half[j]])
				maxc = rem_half[j];
			else if (subscore[t][maxc] == subscore[t][rem_half[j]])
			{
				if (set_age[t][maxc] < set_age[t][rem_half[j]])
					maxc = rem_half[j];
			}
		}
	}
	sad = maxc;
}

void remove(int t, int c)
{
	cs[t][c].is_in_c = 0;
	cs[t][c].score = -cs[t][c].score;
	cs[t][c].config = 0;
	set_age[t][c] = 0.9*set_age[t][c];
	int i, j, h, l;
	for (h = 0; h < nrow[t][c]; h++)
	{
		i = col[t][c][h];
		covered[t][i] --;
		if (covered[t][i] == 0)
		{
			for (l = 0; l < ncol[t][i]; l++)
			{
				j = row[t][i][l];
				if (j == c)
					continue;
				cs[t][j].config = 1;
				cs[t][j].score += rs[t][i].weight;
				subscore[t][j] = subscore[t][j] - rs[t][i].weight;
			}
			index_uncover[t][i] = 1;
			uncover_num++;
		}
		else if (covered[t][i] == 1)
		{
			for (l = 0; l < ncol[t][i]; l++)
			{
				j = row[t][i][l];
				if (j == c)
					continue;
				cs[t][j].config = 1;
				if (cs[t][j].is_in_c == 0) {
					subscore[t][j] += rs[t][i].weight;
				}
				else {
					subscore[t][j] -= rs[t][i].weight;
					cs[t][j].score -= rs[t][i].weight;
				}
			}
		}
		else if (covered[t][i] == 2)
		{
			for (l = 0; l < ncol[t][i]; l++)
			{
				j = row[t][i][l];
				if (j == c)
					continue;
				cs[t][j].config = 1;
				if (cs[t][j].is_in_c == 1) {
					subscore[t][j] += rs[t][i].weight;
				}
			}
		}
		else if (covered[t][i] > 2) {
			for (l = 0; l < ncol[t][i]; l++)
			{
				j = row[t][i][l];
				if (j == c)
					continue;
				cs[t][j].config = 1;
			}
		}
	}
}
void find_ran_vout(int t)
{
	int maxc = 0, varr = 0;
	int maxx = INT_MIN;
	int count1 = 0, ram = 0;
	ram = rand() % uncover_num;
	for (int j = 0; j < m; j++)
	{
		if (index_uncover[t][j] == -1)
			continue;
		if (count1 == ram)
		{
			varr = j;
			break;
		}
		count1++;
	}
	//8.19 randstep
	int therand = ncol[t][varr];
	randstep = rand() % therand;
	//
	int jj = 0;
	int jpj = 0;
	for (int j = randstep; j < ncol[t][varr]; j++)
	{
		jj = row[t][varr][j];
		if (cs[t][jj].is_in_c == 1)
			continue;
		jpj = 1;
		if (maxx < cs[t][jj].score)
		{
			maxc = jj;
			maxx = cs[t][jj].score;
		}
		else if (maxx == cs[t][jj].score)
		{
			if (subscore[t][maxc] < subscore[t][jj])
				maxc = jj;
			else if (subscore[t][maxc] == subscore[t][jj]) {
				if (set_age[t][maxc] < set_age[t][jj])
					maxc = jj;
			}
		}
	}
	for (int j = 0; j < randstep; j++)
	{
		jj = row[t][varr][j];
		if (cs[t][jj].is_in_c == 1)
			continue;
		jpj = 1;
		if (maxx < cs[t][jj].score)
		{
			maxc = jj;
			maxx = cs[t][jj].score;
		}
		else if (maxx == cs[t][jj].score)
		{
			if (subscore[t][maxc] < subscore[t][jj])
				maxc = jj;
			else if (subscore[t][maxc] == subscore[t][jj]) {
				if (set_age[t][maxc] < set_age[t][jj])
					maxc = jj;
			}
		}
	}
	happy = maxc;
}




void find_random_in_vout(int t)
{
	int maxc = 0, varr = 0;
	int maxx = INT_MIN;
	int count1 = 0, ram = 0;
	ram = rand() % uncover_num;
	for (int j = 0; j < m; j++)
	{
		if (index_uncover[t][j] == -1)
			continue;
		if (count1 == ram)
		{
			varr = j;
			break;
		}
		count1++;
	}
	int therand = ncol[t][varr];
	randstep = rand() % therand;
	//
	int jj = 0;
	int jpj = 0;
	for (int j = randstep; j < ncol[t][varr]; j++)
	{
		jj = row[t][varr][j];
		if (cs[t][jj].is_in_c == 1 || cs[t][jj].config == 0)
			continue;
		jpj = 1;
		if (maxx < cs[t][jj].score)
		{
			maxc = jj;
			maxx = cs[t][jj].score;
		}
		else if (maxx == cs[t][jj].score)
		{
			if (subscore[t][maxc] < subscore[t][jj])
				maxc = jj;
			else if (subscore[t][maxc] == subscore[t][jj]) {
				if (set_age[t][maxc] < set_age[t][jj])
					maxc = jj;
			}
		}
	}
	for (int j = 0; j < randstep; j++)
	{
		jj = row[t][varr][j];
		if (cs[t][jj].is_in_c == 1 || cs[t][jj].config == 0)
			continue;
		jpj = 1;
		if (maxx < cs[t][jj].score)
		{
			maxc = jj;
			maxx = cs[t][jj].score;
		}
		else if (maxx == cs[t][jj].score)
		{
			//subscore
			if (subscore[t][maxc] < subscore[t][jj])
				maxc = jj;
			else if (subscore[t][maxc] == subscore[t][jj]) {
				if (set_age[t][maxc] < set_age[t][jj])
					maxc = jj;
			}
		}
	}
	happy = maxc;
	int pjp = 0;
	if (jpj == 0) {
		if (time_test.NeedToStop())	return;

		pjp++;
		if (pjp > max_k)
		{
			find_ran_vout(t);
		}
		else
			find_random_in_vout(t);
	}
}

void add(int t, int c)
{
	cs[t][c].is_in_c = 1;
	cs[t][c].score = -cs[t][c].score;
	set_age[t][c] = 0.9*set_age[t][c];
	tabu = c;
	int i, j, h, l;
	for (h = 0; h < nrow[t][c]; h++)
	{
		i = col[t][c][h];
		covered[t][i] ++;

		if (covered[t][i] == 1)
		{
			for (l = 0; l < ncol[t][i]; l++)
			{
				j = row[t][i][l];
				if (j == c)
					continue;
				cs[t][j].config = 1;//config = 1;
				cs[t][j].score -= rs[t][i].weight;
				subscore[t][j] += rs[t][i].weight;
			}
			index_uncover[t][i] = -1;
			uncover_num--;
		}
		else if (covered[t][i] == 2)
		{
			for (l = 0; l < ncol[t][i]; l++)
			{
				j = row[t][i][l];
				if (j == c)
					continue;
				cs[t][j].config = 1;
				if (cs[t][j].is_in_c == 0) {
					subscore[t][j] -= rs[t][i].weight;
				}
				else {
					cs[t][j].score += rs[t][i].weight;
					subscore[t][j] += rs[t][i].weight;
				}
			}
		}
		else if (covered[t][i] == 3)
		{
			for (l = 0; l < ncol[t][i]; l++)
			{
				j = row[t][i][l];
				if (j == c)
					continue;
				cs[t][j].config = 1;
				if (cs[t][j].is_in_c == 1) {
					subscore[t][j] -= rs[t][i].weight;
				}
			}
		}
		else if (covered[t][i] > 3)
		{
			for (l = 0; l < ncol[t][i]; l++)
			{
				j = row[t][i][l];
				if (j == c)
					continue;
				cs[t][j].config = 1;
			}
		}
	}
}

void reload_weight_score(int t)
{
	
	int fpx, fpx1, fpx2;
	int r = 1000;
	double p = 0.3;
	int init_weight = 0;
	for (fpx = 0; fpx < m; fpx++)
	{
		if (index_uncover[t][fpx] == -1 || fpx == sad)
			continue;
		if (rs[t][fpx].weight > r)
		{
			init_weight = rs[t][fpx].weight * p;
			for (fpx1 = 0; fpx1 < ncol[t][fpx]; fpx1++)
			{
				fpx2 = row[t][fpx][fpx1];
				cs[t][fpx2].score = cs[t][fpx2].score - rs[t][fpx].weight + init_weight;
			}
			rs[t][fpx].weight = init_weight;
		}
		else
		{
			rs[t][fpx].weight = rs[t][fpx].weight + 1;
			for (fpx1 = 0; fpx1 < ncol[t][fpx]; fpx1++)
			{
				fpx2 = row[t][fpx][fpx1];
				cs[t][fpx2].score = cs[t][fpx2].score + 1;

			}
		}
	}
}


void reload_age(int t)
{
	for (int traversal = 0; traversal < n; traversal++)
	{
		if (traversal != happy && traversal != sad)
		{
			set_age[t][traversal] ++;
		}
	}
}


int del_remtemp()
{
	int i;
	for (i = 0; i < max_k; i++)
	{
		if (rem_temp[i] == sad)
			break;
	}
	return i;
}

void CCA_local_search(int t)
{

	int remtemp = 0;

	find_worst_in_c(t);
	remove(t, sad);
	remtemp = del_remtemp();
	find_random_in_vout(t);
	add(t, happy);

	rem_temp[remtemp] = happy;
	reload_age(t);

}

void CCA_local_searchw(int t)
{
	int remtemp = 0;
	find_worst_in_c(t);
	remove(t, sad);
	remtemp = del_remtemp();
	find_random_in_vout(t);
	add(t, happy);
	rem_temp[remtemp] = happy;
	reload_weight_score(t);
	reload_age(t);

}


void init_aco()
{
	int t, i, j;
	for (t = 0; t < ants; t++) {
		for (j = 0; j < n; j++)
			uncover[t][j] = nrow[t][j];
	}
	for (j = 0; j < n; j++)
	{
		tau[j] = 0.5;
		eta[j] = 1.0*uncover[0][j];
	}

	for (j = 0; j < n; j++)
	{
		te_beta[j] = tau[j];
		for (i = 0; i < beta; i++)
			te_beta[j] = te_beta[j] * eta[j];
	}
}


void init_eta_tebeta()
{
	int i, j;
	for (i = 0; i < MAXA; i++)
	{
		for (j = 0; j < n; j++)
			uncover[i][j] = nrow[i][j];
	}

	for (j = 0; j < n; j++)
	{
		eta[j] = 1.0*nrow[0][j];
		te_beta[j] = tau[j];
		for (i = 0; i < beta; i++)
			te_beta[j] = te_beta[j] * eta[j];
	}

}

void add_aoc(int t, int c)
{
;
	int i, j, h, l;
	te_beta[c] = 0.0;
	eta[c] = 0.0;
	cs[t][c].is_in_c = 1;
	cs[t][c].score = -cs[t][c].score;
	for (h = 0; h < nrow[t][c]; h++)
	{
		i = col[t][c][h];
		covered[t][i] ++;
		if (covered[t][i] == 1)
		{
			for (l = 0; l < ncol[t][i]; l++)
			{
				j = row[t][i][l];
				if (j == c)
					continue;
				cs[t][j].config = 1;
				cs[t][j].score -= rs[t][i].weight;
				subscore[t][j] += rs[t][i].weight;
				if (eta[j] != 0.) {
					uncover[t][j] --;
					eta[j] = uncover[t][j] / sqrt(m);
				}
				te_beta[j] = tau[j];
				for (int t = 0; t < beta; t++)
					te_beta[j] = te_beta[j] * eta[j];
			}
			index_uncover[t][i] = -1;
			uncover_num--;
		}
		if (covered[t][i] == 2)
		{
			for (l = 0; l < ncol[t][i]; l++)
			{
				j = row[t][i][l];
				if (j == c)
					continue;
				cs[t][j].config = 1;
				if (cs[t][j].is_in_c == 0) {
					subscore[t][j] -= rs[t][i].weight;
				}
				else {
					cs[t][j].score += rs[t][i].weight;
					subscore[t][j] += rs[t][i].weight;
				}
			}
		}
		if (covered[t][i] == 3)
		{
			for (l = 0; l < ncol[t][i]; l++)
			{
				j = row[t][i][l];
				if (j == c)
					continue;
				cs[t][j].config = 1;
				if (cs[t][j].is_in_c == 1) {
					subscore[t][j] -= rs[t][i].weight;
				}
			}
		}
		if (covered[t][i] > 3)
		{
			for (l = 0; l < ncol[t][i]; l++)
			{
				j = row[t][i][l];
				if (j == c)
					continue;
				cs[t][j].config = 1;
			}
		}
	}
}

void reload_alltau()
{
	for (int i = 0; i < n; i++)
	{
		tau[i] = tau[i] * omiga;
	}
}

void reload_alltaupro(int c)
{
	if (tau[c] < taumin)
		tau[c] = taumin;
	if (tau[c] > taumax)
		tau[c] = taumax;
}

void reload_bt(int c)
{
	tau[c] = tau[c] + delta_tau / 3;
}


void reload_bs(int c)
{
	tau[c] = tau[c] + (delta_tau * 2) / 3;//更新该集合的信息素tau
}





void maxmin()
{
	taumax = best_value / (m * omiga);
	taumin = 1 / (sqrt(m) * sqrt(sqrt(n)));
}


void add_rand(int t, int randi)
{
	int i, mid, randii;
	int pre_k = max_k;
	for (i = 0; i < randi; i++)
	{
		randii = rand() % pre_k;
		add_aoc(t, Remember[randii]);
		pre_k--;
		mid = Remember[randii];
		Remember[randii] = Remember[pre_k];
		Remember[pre_k] = mid;
	}
}
void ANT_1()
{
	int t = 0;
	init_eta_tebeta();
	ant_empty(t);
	for (int i = 0; i < max_k; i++)
	{
		add_aoc(t, Remember[i]);
	}
	ant_search(t);
}

void ANT_2()
{
	int t = 1;
	init_eta_tebeta();
	ant_empty(t);
	ant_build(t, max_k);
	ant_search(t);
}

void ANT_3()
{
	int t = 2;
	init_eta_tebeta();
	ant_empty(t);
	//0%--50%
	randant = rand() % valuek;
	add_rand(t, randant);
	int itor = max_k - randant;
	ant_build(t, itor);
	ant_search(t);
}
void ANT_4()
{
	int t = 3;
	init_eta_tebeta();
	ant_empty(t);
	//50%--
	randant = rand() % valuek;
	int A = max_k / 2 + randant;
	add_rand(t, A);
	int itor = max_k - A;
	ant_build(t, itor);
	ant_search(t);
}
void ANT_5()
{
	int t = 4;
	init_eta_tebeta();
	ant_empty(t);
	for (int i = 0; i < max_k; i++)
	{
		add_aoc(t, Remember[i]);
	}
	ant_searchw(t);
}
void ANT_6()
{
	int t = 5;
	init_eta_tebeta();
	ant_empty(t);
	ant_build(t, max_k);
	ant_searchw(t);
}

void ANT_7()
{
	int t = 6;
	init_eta_tebeta();
	ant_empty(t);
	//0%-50%
	randant = rand() % valuek;
	add_rand(t, randant);
	int itor = max_k - randant;
	ant_build(t, itor);
	ant_searchw(t);
}
void ANT_8()
{
	int t = 7;
	init_eta_tebeta();
	ant_empty(t);
	//50%--
	randant = rand() % valuek;
	int A = max_k / 2 + randant;
	add_rand(t, A);
	int itor = max_k - A;
	ant_build(t, itor);
	ant_searchw(t);
}

void ANT_greedy()
{
	int t = 10;
	init_eta_tebeta();
	ant_empty_greedy(t);
	for (int i = 0; i < max_k; i++)
	{
		add_aoc_greedy(t, Remember[i]);
	}
	ant_search_greedy(t);
}


void ant_empty_greedy(int t)
{
	for (int i = 0; i < n; i++)
	{
		if (cs[t][i].is_in_c == 1)
		{
			remove_greedy(t, i);
		}
	}
	uncover_num = m;
}

void remove_greedy(int t, int c)
{
	cs[t][c].is_in_c = 0;
	cs[t][c].score = -cs[t][c].score;
	int i, j, h, l;
	for(h = 0; h < nrow[t][c]; h++)
	{
		i = col[t][c][h];
		covered[t][i] --;
		if (covered[t][i] == 0)
		{
			for (l = 0; l < ncol[t][i]; l++)
			{
				j = row[t][i][l];
				if (j == c)
					continue;
				cs[t][j].score += rs[t][i].weight;
			}
			index_uncover[t][i] = 1;
			uncover_num++;
		}
		else if (covered[t][i] == 1)
		{
			for (l = 0; l < ncol[t][i]; l++)
			{
				j = row[t][i][l];
				if (j == c)
					continue;
				if (cs[t][j].is_in_c != 0) {
					cs[t][j].score -= rs[t][i].weight;
				}
			}
		}
	}
}

void add_aoc_greedy(int t, int c)
{
	int i, j, h, l;
	te_beta[c] = 0.0;
	eta[c] = 0.0;
	cs[t][c].is_in_c = 1;
	cs[t][c].score = -cs[t][c].score;
	for (h = 0; h < nrow[t][c]; h++)
	{
		i = col[t][c][h];
		covered[t][i] ++;
		if (covered[t][i] == 1)
		{
			for (l = 0; l < ncol[t][i]; l++)
			{
				j = row[t][i][l];
				if (j == c)
					continue;
				cs[t][j].score -= rs[t][i].weight;
				if (eta[j] != 0.) {
					uncover[t][j] --;
					eta[j] = uncover[t][j] / sqrt(m);
				}
				te_beta[j] = tau[j];
				for (int t = 0; t < beta; t++)
					te_beta[j] = te_beta[j] * eta[j];
			}
			index_uncover[t][i] = -1;
			uncover_num--;
		}
		if (covered[t][i] == 2)
		{
			for (l = 0; l < ncol[t][i]; l++)
			{
				j = row[t][i][l];
				if (j == c)
					continue;
				if (cs[t][j].is_in_c != 0) {
					cs[t][j].score += rs[t][i].weight;
				}
			}
		}
	}
}

void add_greedy(int t,int c)
{
	cs[t][c].is_in_c = 1;
	cs[t][c].score = -cs[t][c].score;
	//tabu start
	tabu = c;
	//tabu end
	int i, j, h, l;
	for (h = 0; h < nrow[t][c]; h++)
	{
		i = col[t][c][h];
		covered[t][i] ++;
		if (covered[t][i] == 1)
		{
			for (l = 0; l < ncol[t][i]; l++)
			{
				j = row[t][i][l];
				if (j == c)
					continue;
				cs[t][j].score -= rs[t][i].weight;
			}
			index_uncover[t][i] = -1;
			uncover_num--;
		}
		else if (covered[t][i] == 2)
		{
			for (l = 0; l < ncol[t][i]; l++)
			{
				j = row[t][i][l];
				if (j == c)
					continue;
				if (cs[t][j].is_in_c != 0) {
					cs[t][j].score += rs[t][i].weight;
				}
			}
		}
	}	
}

void ant_search_greedy(int t)
{
	cntants = 0;
	tpp = 0;
	int r = 0;
	tpp = 0;
	for (r = 0; r < n; r++)
	{
		if (cs[t][r].is_in_c == 1)
		{
			rem_temp[tpp] = r;
			tpp++;
		}
	}
	while (1)
	{
		greedy_ls(t);

		if (m - uncover_num > best_value)
		{
			tpp = 0;
			cntants = 0;
			for (r = 0; r < n; r++)
			{
				if (cs[t][r].is_in_c == 1)
				{
					Remember[tpp] = r;
					tpp++;
				}
			}
			update_best_sol();
		}
		else
		{
			cntants++;
			if (cntants > greedy_limit)
			{
				uncover_num = uncover_temp;
				break;
			}
		}
		if (time_test.NeedToStop())	return;
	}
}

void greedy_ls(int t)
{
	int remtemp = 0;
	greedy_worst(t);
	remove_greedy(t, sad);
	remtemp = del_remtemp();
	greedy_in_vout(t);
	add(t, happy);
	rem_temp[remtemp] = happy;
}


void greedy_worst(int t)
{
	int maxc = 0, maxx = INT_MIN;
	int k_temp = max_k;
	int tempk = 0;
	for (int j = 0; j < max_k; j++)
	{
		randnum = rand() % k_temp;
		rem_half[j] = rem_temp[randnum];
		k_temp--;
		tempk = rem_temp[randnum];
		rem_temp[randnum] = rem_temp[k_temp];
		rem_temp[k_temp] = tempk;
	}
	randstep = rand() % k_half;
	for (int j = randstep; j < k_half; j++)
	{
		if (rem_half[j] == tabu)
			continue;
		if (maxx < cs[t][rem_half[j]].score)
		{
			maxc = rem_half[j];
			maxx = cs[t][rem_half[j]].score;
		}
	}
	for (int j = 0; j < randstep; j++)
	{

		if (rem_half[j] == tabu)
			continue;

		if (maxx < cs[t][rem_half[j]].score)
		{
			maxc = rem_half[j];
			maxx = cs[t][rem_half[j]].score;
		}
	}
	sad = maxc;
}

void greedy_in_vout(int t)
{
	int maxc = 0, varr = 0;
	int maxx = INT_MIN;
	int count1 = 0, ram = 0;
	ram = rand() % uncover_num;
	for (int j = 0; j < m; j++)
	{
		if (index_uncover[t][j] == -1)
			continue;
		if (count1 == ram)
		{
			varr = j;
			break;
		}
		count1++;
	}
	int therand = ncol[t][varr];
	randstep = rand() % therand;
	int jj = 0;
	int jpj = 0;
	for (int j = randstep; j < ncol[t][varr]; j++)
	{
		jj = row[t][varr][j];
		if (cs[t][jj].is_in_c == 1)
			continue;
		jpj = 1;
		if (maxx < cs[t][jj].score)
		{
			maxc = jj;
			maxx = cs[t][jj].score;
		}
	}
	for (int j = 0; j < randstep; j++)
	{
		jj = row[t][varr][j];
		if (cs[t][jj].is_in_c == 1)
			continue;
		jpj = 1;
		if (maxx < cs[t][jj].score)
		{
			maxc = jj;
			maxx = cs[t][jj].score;
		}
	}
	happy = maxc;
}

void ant_empty(int t)
{
	for (int i = 0; i < n; i++)
	{
		if (cs[t][i].is_in_c == 1)
		{
			remove(t, i);
		}
	}
	uncover_num = m;
}


void ant_build(int t, int itor)
{
	int antk;
	for (antk = 0; antk < itor; antk++)
	{
		int j = 0;
		double pper = 0.;
		int Irand = rand() % INT_MAX;
		double Frand = Irand / (INT_MAX * 1.0);
		double temp_sum = 0.;
		for (j = 0; j < n; j++)
		{
			if (cs[t][j].is_in_c == 1)
				continue;
			temp_sum = temp_sum + te_beta[j];
		}
		for (j = 0; j < n; j++)
		{
			if (cs[t][j].is_in_c == 0)
			{
				pper += te_beta[j] / temp_sum;
			}
			if (pper > Frand)
				break;
		}
		add_aoc(t, j);
	}
}

void ant_search(int t)
{
	cntants = 0;
	tpp = 0;
	int r = 0;
	for (r = 0; r < n; r++)
	{
		if (cs[t][r].is_in_c == 1)
		{
			rem_temp[tpp] = r;
			tpp++;
		}
	}
	while (1)
	{
		CCA_local_search(t);

		if (m - uncover_num > best_value)
		{
			tpp = 0;
			cntants = 0;
			for (r = 0; r < n; r++)
			{
				if (cs[t][r].is_in_c == 1)
				{
					Remember[tpp] = r;
					tpp++;
				}
			}
			update_best_sol();
		}
		else
		{
			cntants++;
			if (cntants > cca_limit)
			{
				uncover_num = uncover_temp;
				break;
			}
		}
		if (time_test.NeedToStop())	return;
	}
}

void ant_searchw(int t)
{
	cntants = 0;
	tpp = 0;
	int r = 0;
	tpp = 0;
	for (r = 0; r < n; r++)
	{
		if (cs[t][r].is_in_c == 1)
		{
			rem_temp[tpp] = r;
			tpp++;
		}
	}
	while (1)
	{
		CCA_local_searchw(t);
		if (m - uncover_num > best_value)
		{
			tpp = 0;
			cntants = 0;
			for (r = 0; r < n; r++)
			{
				if (cs[t][r].is_in_c == 1)
				{
					Remember[tpp] = r;
					tpp++;
				}
			}
			update_best_sol();
		}
		else
		{
			cntants++;
			if (cntants > cca_limit)
			{
				uncover_num = uncover_temp;
				break;
			}
		}
		if (time_test.NeedToStop())	return;
	}
}

void init_search(int t)
{
	cntants = 0;
	tpp = 0;
	int r = 0;
	tpp = 0;
	for (r = 0; r < n; r++)
	{
		if (cs[t][r].is_in_c == 1)
		{
			rem_temp[tpp] = r;
			tpp++;
		}
	}
	while (1)
	{
		CCA_local_search(t);

		if (m - uncover_num > best_value)
		{
			tpp = 0;
			cntants = 0;
			for (r = 0; r < n; r++)
			{
				if (cs[t][r].is_in_c == 1)
				{
					Remember[tpp] = r;
					tpp++;
				}
			}
			update_best_sol();
		}
		else
		{
			cntants++;
			if (cntants > 100)
			{
				uncover_num = uncover_temp;
				break;
			}
		}
		if (time_test.NeedToStop())	return;
	}
}

void init_searchw(int t)
{
	cntants = 0;
	tpp = 0;
	int r = 0;
	tpp = 0;
	for (r = 0; r < n; r++)
	{
		if (cs[t][r].is_in_c == 1)
		{
			rem_temp[tpp] = r;
			tpp++;
		}
	}
	while (1)
	{
		CCA_local_searchw(t);
		if (m - uncover_num > best_value)
		{
			tpp = 0;
			cntants = 0;
			for (r = 0; r < n; r++)
			{
				if (cs[t][r].is_in_c == 1)
				{
					Remember[tpp] = r;
					tpp++;
				}
			}
			update_best_sol();
		}
		else
		{
			cntants++;
			if (cntants > 100)
			{
				uncover_num = uncover_temp;
				break;
			}
		}
		if (time_test.NeedToStop())	return;
	}
}
int adc = 0;
void MMAS_pro_density()
{
	valuek = max_k / 5;
	init_eta_tebeta();
	if (tag_init == 1)
	{
		int j = 0;
		tag_init = 0;
		uncover_num = m;
		int tempi = 8;
		FourlevelGreedy(tempi);
		for (j = 0; j < max_k; j++)
		{
			remember[j] = Remember[j];
		}
		best_temp = best_value;
		best_value = max_k;
		ANT_1();
		if (finish_time > time_limit)
			return;
		tempi = 9;
		uncover_num = m;
		FourlevelGreedyw(tempi);
		ANT_5();
		if (finish_time > time_limit)
			return;

		if (best_value < best_temp)
		{
			best_value = best_temp;
			best_temp = best_value;
		}

		for (j = 0; j < max_k; j++)
		{
			rtemp[j] = Remember[j];
			Remember[j] = remember[j];
		}
		ANT_2();
		if (finish_time > time_limit)
			return;
		ANT_4();
		if (finish_time > time_limit)
			return;
		for (j = 0; j < max_k; j++)
		{
			remember[j] = Remember[j];
			Remember[j] = rtemp[j];
		}
		int lintemp = best_value;
		best_value = best_temp;
		best_temp = lintemp;
		ANT_6();
		if (finish_time > time_limit)
			return;
		ANT_8();
		if (finish_time > time_limit)
			return;

		if (best_value < best_temp)
		{
			best_value = best_temp;
			best_temp = best_value;
			for (j = 0; j < max_k; j++)
			{
				Remember[j] = rtemp[j];
			}
		}
	}
	else if (tag_init == 0)
	{
		ANT_1();
		if (finish_time > time_limit)
			return;
		ANT_5();
		if (finish_time > time_limit)
			return;
		ANT_2();
		if (finish_time > time_limit)
			return;
		ANT_4();
		if (finish_time > time_limit)
			return;
		ANT_6();
		if (finish_time > time_limit)
			return;
		ANT_8();
		if (finish_time > time_limit)
			return;
		ANT_greedy();
		if (finish_time > time_limit)
			return;
		
		if(greedy_limit > 10000 && smooth%2 == 1){
			smooth ++;
			greedy_limit = greedy_limit - 30000;
		}
		else
			greedy_limit = 1000;
		if(cca_limit < 10000)
			cca_limit = cca_limit*10;
		else if(cca_limit < 50000)
			cca_limit = cca_limit + 20000;
		//>100,/5;<10000,*10;<50000,+10000
		adc++;
	}
    maxmin();
	reload_alltau();
	for (int j = 0; j < max_k; j++)
	{
		reload_bs(Remember[j]);
		reload_alltaupro(j);
	}
}

void MMAS_pro()
{
	valuek = max_k / 5;
	init_eta_tebeta();
	if (tag_init == 1)
	{
		int j = 0;
		tag_init = 0;
		uncover_num = m;
		int tempi = 8;
		FourlevelGreedy(tempi);
		for (j = 0; j < max_k; j++)
		{
			remember[j] = Remember[j];
		}
		best_temp = best_value;
		best_value = max_k;
		ANT_1();
		if (finish_time > time_limit)
			return;
		tempi = 9;
		uncover_num = m;
		FourlevelGreedyw(tempi);
		ANT_5();
		if (finish_time > time_limit)
			return;

		if (best_value < best_temp)
		{
			best_value = best_temp;
			best_temp = best_value;
		}

		for (j = 0; j < max_k; j++)
		{
			rtemp[j] = Remember[j];
			Remember[j] = remember[j];
		}
		ANT_2();
		if (finish_time > time_limit)
			return;
		ANT_3();
		if (finish_time > time_limit)
			return;
		ANT_4();
		if (finish_time > time_limit)
			return;

		for (j = 0; j < max_k; j++)
		{
			remember[j] = Remember[j];
			Remember[j] = rtemp[j];
		}

		int lintemp = best_value;
		best_value = best_temp;
		best_temp = lintemp;

		ANT_6();
		if (finish_time > time_limit)
			return;

		ANT_7();
		if (finish_time > time_limit)
			return;

		ANT_8();
		if (finish_time > time_limit)
			return;

		if (best_value < best_temp)
		{
			best_value = best_temp;
			best_temp = best_value;
			for (j = 0; j < max_k; j++)
			{
				Remember[j] = rtemp[j];
			}
		}
	}
	else if (tag_init == 0)
	{
		ANT_1();
		if (finish_time > time_limit)
			return;
		ANT_5();
		if (finish_time > time_limit)
			return;
		ANT_2();
		if (finish_time > time_limit)
			return;
		ANT_3();
		if (finish_time > time_limit)
			return;
		ANT_4();
		if (finish_time > time_limit)
			return;
		ANT_6();
		if (finish_time > time_limit)
			return;
		ANT_7();
		if (finish_time > time_limit)
			return;
		ANT_8();
		if (finish_time > time_limit)
			return;
	}
        maxmin();
	reload_alltau();

	for (int j = 0; j < max_k; j++)
	{
		reload_bs(Remember[j]);
		reload_alltaupro(j);
	}
}


void NuAnts_density()
{
	best_value = max_k;
	init();
	init_aco();
	best_temp = best_value;
	maxmin();
	for (int mnb = 0; mnb < n; mnb++)
		tau[mnb] = (taumax + taumin) / 2;
	while (1)
	{
		MMAS_pro_density();
		if (time_test.NeedToStop())	break;
	}
	return;
}


void NuAnts()
{
	best_value = max_k;
	init();
	init_aco();
	best_temp = best_value;
	maxmin();
	for (int mnb = 0; mnb < n; mnb++)
		tau[mnb] = (taumax + taumin) / 2;
	while (1)
	{
		MMAS_pro();
		if (time_test.NeedToStop())	break;
	}
	return;
}


void update_best_sol()
{
	int i, j;
	if (m - uncover_num > best_value)
	{
		real_time = time_test.GetTime();
		best_value = m - uncover_num;
	}
}

int main(int argc, char *argv[])
{
	build_instance(argv[1]);
	max_k = atoi(argv[2]);
	seed = atoi(argv[3]);
	time_limit = atof(argv[4]);
	printf("MMAS %s  K=%d   seed=%d   t=%fs\n", argv[1], max_k, seed, time_limit);
	delta_tau = 1 / sqrt(m);
	srand(seed);
	k_half = (max_k * 2) / 3;
	if( density == 1)
	{
		time_test.Init();
		NuAnts_density();
	}
	else{
		time_test.Init();
		NuAnts();
	}
	freopen("mmas.txt", "a", stdout);
	printf("Instance: %s\t\t K:%d\t\t seed:%d\t\t Result:%d\t\t Time:%.2f\n", argv[1], max_k, seed, best_value, real_time);
	free_all();
	return 0;
}

