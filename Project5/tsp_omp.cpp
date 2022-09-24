#include<bits/stdc++.h>
#include <ctime>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include<omp.h>

//Use this for random number on every run
//#include <time.h> //srand(time(0));

using namespace std;

const int MAXNCITIES = 16; //Maximum Number of City
const int INFINITE = 1e9;
const int MIN_PATH_WEIGHT = 1;
const int MAX_PATH_WEIGHT = 10;

long long factorial[MAXNCITIES + 1];

int randomNum(int min, int max)
{
	return min + rand() % (max - min + 1);
}

//No use (Exp: 10!)
void precompute_factorial()
{
	factorial[0] = 1;

	for (int i = 1; i <= MAXNCITIES; i++)
	{
		factorial[i] = i * factorial[i - 1];
	}
}

//Assign Matrix Values to a Symmetric Matrix
void assign_edge_weights(vector<vector<int>>& matrix)
{
	int n = matrix.size();
	//Use this for random number on every run
	//srand(time(0));

	//Assign Matrix Values to a Symmetric Matrix
	for (int i = 0; i < n; i++)
	{
		for (int j = i + 1; j < n; j++)
		{
			//Assign Matrix
			matrix[i][j] = randomNum(MIN_PATH_WEIGHT, MAX_PATH_WEIGHT);
			matrix[j][i] = matrix[i][j];
		}
		matrix[i][i] = 0;
	}
}

//Display Matrix
void display_matrix(vector<vector<int>>& matrix)
{
	int N = matrix.size();

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			printf("%4d", matrix[i][j]);
		}
		cout << endl;
	}
}

// Find the total current path cost
int find_path_cost(vector<vector<int>>& matrix, vector<int>& arr)
{
	//If 10 cities are inputted then arr is 11, so loop 10 times
	int current_cost = 0;

	for (int i = 1; i < (int)arr.size(); i++)
	{
		current_cost += matrix[arr[i]][arr[i - 1]];
	}
	return current_cost;
}

// Here Can Be Implement In Parallel Algorithm
vector<int> nth_permutation(vector<int>& arr, long long n)
{
	int N = arr.size();

	assert(n <= factorial[N]);

	sort(arr.begin(), arr.end());

	set<int>st;
	for (int x : arr)st.insert(x);

	vector<int>ans;


	for (int i = 0; i < N; i++)
	{
		int cn = 1;
		long long cval = factorial[N - 1 - i];
		while (cval < n)
		{
			cn++;
			cval = (long long)cn * cval;
			cval = (long long)cval / (cn - 1);
		}

		long long pval = cval * (cn - 1) / cn;
		n -= pval;

		auto it = st.begin();

		for (int i = 0; i < cn - 1; i++)it++;
		ans.push_back(*it);
		st.erase(it);
	}

	return ans;
}

vector<int> tsp_omp(vector<vector<int>>& matrix)
{
	int n = matrix.size();

	int optimal_value = INFINITE; //it is an algorithm for infinity (optimal_cost = infinity)

	vector<int>route;
	vector<int>nodes;

	long int k = factorial[n - 1];

#pragma omp parallel private(nodes) shared(route,optimal_value)
	{

		for (int i = 1; i < n; i++)nodes.push_back(i); //If 10 cities then push 1-9, Exp: 1, 2, 3, ..., 9

		////Check the initial nodes
		//for (int i = 0; i < nodes.size(); i++) {
		//	printf("%d", nodes[i]);
		//}

		int num = omp_get_num_threads();
		int id = omp_get_thread_num();
		long int iter_per_thread = k / num;
		int extra = k % num;

		if (id < extra) {
			nodes = nth_permutation(nodes, (id) * (iter_per_thread + 1));
			iter_per_thread = iter_per_thread + 1;
		}
		else nodes = nth_permutation(nodes, (id)*iter_per_thread + extra);

		int i = 0;
		do
		{
			//Assign nodes to temporary
			vector<int>tempNodes = nodes;

			//Push last nodes is 0
			tempNodes.push_back(0);

			//Assign first nodes is 0
			tempNodes.insert(tempNodes.begin(), 0);

			/*for (int i = 0; i < temp.size(); i++) {
			printf("%d\n", temp[i]);
			}
			printf("\n");
			display_matrix(matrix);*/

			//CurrentPCost is to get the current total path cost of the random matrix
			int currentPCost = find_path_cost(matrix, tempNodes);
#pragma omp critical
			{
				if (currentPCost < optimal_value) {
					optimal_value = currentPCost;
					route = tempNodes;
				}
			}
			i++;

			next_permutation(nodes.begin(), nodes.end());
		} while (i < iter_per_thread);
	}

	return route;
}

int main(int argc, char* argv[])
{
	if (argc != 3) {
		std::cout << "Please Input A Proper Argument (exe cities threads) (Exp: Project5.exe 10 4)" << std::endl;
		return -1;
	}

	int N = stoi(argv[1]); //stoi is for convert string to integer
	int threads = stoi(argv[2]); //get user number of thread

	omp_set_num_threads(threads);

	precompute_factorial();

	//Run a 2D Matrix, N is Row, vector<int>(N, 0) is column, Start a Matrix with 0 values
	vector<vector<int>>matrix(N, vector<int>(N, 0));
	assign_edge_weights(matrix);

	// Display the path weight matrix
	printf("Path Weight Matrix: \n");
	display_matrix(matrix);
	cout << endl;

	auto start = std::chrono::high_resolution_clock::now();		// start time

	vector<int>route = tsp_omp(matrix);

	auto finish = std::chrono::high_resolution_clock::now();	// end time

	// Display the minimum cost path
	printf("Minimum Cost Path for %d cities: \n", N);
	for (auto x : route)cout << x << " "; //no need put auto range-declaration
	cout << endl;
	cout << endl;

	// Display the minimum path cost
	printf("Minimum Total Path Cost for %d cities: \n", N);
	cout << find_path_cost(matrix, route) << endl;
	cout << endl;

	// Display the run-time
	chrono::duration<double> elapsed = finish - start;
	cout << fixed << setprecision(5) << elapsed.count() << endl;

	printf("Run-Time for %d cities: ", N);
	cout << fixed << setprecision(5) << elapsed.count() << endl;
}
