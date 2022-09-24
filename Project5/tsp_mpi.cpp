#define MASTERPE 0				/* The MASTERPE PE					*/

#include <bits/stdc++.h>
#include <cstdlib>
#include <iostream>
#define NOMINMAX
#include <windows.h>
#include "mpi.h"

//Use this for random number on every run
//#include <time.h> //srand(time(0));

using namespace std;

const int MAXNCITIES = 16;
const int INFINITECOST = 1e9;
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
void assign_edge_weights(int** matrix, int N)
{
	//Use this for random number on every run
	//srand(time(0));

	//Assign Matrix Values to a Symmetric Matrix
	for (int i = 0; i < N; i++) {
		for (int j = i + 1; j < N; j++) {
			matrix[i][j] = randomNum(MIN_PATH_WEIGHT, MAX_PATH_WEIGHT);
			matrix[j][i] = matrix[i][j];
		}
		matrix[i][i] = 0;
	}
}

//Display Matrix
void display_matrix(int** matrix, int N)
{
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
int find_path_cost(int** matrix, vector<int>& arr)
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

	vector<int>reply;

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
		reply.push_back(*it);
		st.erase(it);
	}

	return reply;
}

int main(int argc, char* argv[])
{
	//Set Console OUTPUT Size =========================================================================================================

	HWND console = GetConsoleWindow();
	RECT r;
	GetWindowRect(console, &r); //stores the console's current dimensions
	MoveWindow(console, r.left, r.top, 900, 800, TRUE);

	if (argc != 2) {
		std::cout << "Please Input A Proper Argument (exe cities) (Exp: Project5.exe 10)" << std::endl;
		return -1;
	}

	int N = stoi(argv[1]);

	int			totalPEs;					/* Number of PEs */
	int			currentPE;					/* My PE number  */
	int			errorStatus;					/* Error stat  */
	MPI_Status	stat;

	double	startTime, endTime;					/* timing */

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &totalPEs);
	MPI_Comm_rank(MPI_COMM_WORLD, &currentPE);

	precompute_factorial();

	if (currentPE == MASTERPE) {
		//Main Header OUTPUT ==============================================================================================================
	// 
		HANDLE std_output = GetStdHandle(STD_OUTPUT_HANDLE);
		//Set Color
		SetConsoleTextAttribute(std_output, 11);
		printf("========================================================================================================================");
		printf("\n");
		printf("						 ________ _____ _____");
		printf("\n");
		printf("						| __   __/ ____ | __ \\");
		printf("\n");
		printf("						    | |  | (___ | |__)|");
		printf("\n");
		printf("						    | |   \\___ \\| ___/");
		printf("\n");
		printf("						    | |   ____) | |");
		printf("\n");
		printf("						    |_|  |_____/|_|");
		printf("\n");
		printf("========================================================================================================================");
		printf("\n");
		printf("				Travelling Salesman Problem: Parallel Implementations");
		printf("\n");
		printf("========================================================================================================================");
		printf("\n");
		SetConsoleTextAttribute(std_output, 10);
		printf("				      Message Passing Interface(MPI)");
		printf("\n");
		SetConsoleTextAttribute(std_output, 11);
		printf("========================================================================================================================");
		printf("\n");
		printf("\n");
		printf("\n");
		SetConsoleTextAttribute(std_output, 7);
	}

	int** matrix = new int* [N];
	for (int i = 0; i < N; i++) {
		matrix[i] = new int[N];
	}

	// assigning edge weights in MASTERPE
	if (currentPE == MASTERPE) {
		assign_edge_weights(matrix, N);
	}
	// sending the edge weights to other PEs
	if (currentPE == MASTERPE) {
		for (int pe = 1; pe < totalPEs; pe++) {
			for (int i = 0; i < N; i++) {
				MPI_Send(&matrix[i][0], N, MPI_INT, pe, pe, MPI_COMM_WORLD);
				/*MPI_Send(void* data, int count, MPI_Datatype datatype, int destination, int tag,
				MPI_Comm communicator)*/
				//tag = ID of the message
			}
		}
	}
	else {// receiving the PE from the master PE
		for (int i = 0; i < N; i++) {
			MPI_Recv(&matrix[i][0], N, MPI_INT, MASTERPE, currentPE, MPI_COMM_WORLD, &stat);
			/*MPI_Recv(void* data, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm communicator,
				MPI_Status * stat)*/
		}
	}

	MPI_Barrier(MPI_COMM_WORLD); // to synchronize all processes in a communicator 

	// Display the path weight matrix
	if (currentPE == MASTERPE) {
		printf("Path Weight Matrix: \n");
		display_matrix(matrix, N);
		cout << endl;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	startTime = MPI_Wtime();	// startTime time

	// MPI parallel

	int optimal_cost = INFINITECOST; //it is an algorithm for infinity (optimal_cost = infinity)
	int* reply = new int[N + 1];

	// divide among PEs
	//long takes 64-bit (8-bytes) of space, and the long long takes 128-bits (16-bytes) of space
	long long numberPermutationPEs = factorial[N - 1] / totalPEs; //number of permutations of the PE
	long long remainderOfPermutation = factorial[N - 1] % totalPEs; // the remainder of the permutation

	long long initial_perm_index, end_perm_index;


	vector<int>nodes;
	for (int i = 1; i < N; i++) nodes.push_back(i); //If 10 cities then push 1-9, Exp: 1, 2, 3, ..., 9

	////Check the initial nodes
	//for (int i = 0; i < nodes.size(); i++) {
	//	printf("%d", nodes[i]);
	//}

	//remaind of permutation is 0
	if (remainderOfPermutation == 0) {
		initial_perm_index = (currentPE * numberPermutationPEs) + 1;
		end_perm_index = (currentPE + 1) * numberPermutationPEs;
	}
	else {//remaind of permutation is bigger than the current PE
		if (currentPE < remainderOfPermutation) {
			initial_perm_index = (currentPE * (numberPermutationPEs + 1)) + 1;
			end_perm_index = (currentPE + 1) * (numberPermutationPEs + 1);
		}
		else {//remaind of permutation is less than the current PE
			initial_perm_index = remainderOfPermutation * (numberPermutationPEs + 1) + (currentPE - remainderOfPermutation) * numberPermutationPEs + 1;
			end_perm_index = remainderOfPermutation * (numberPermutationPEs + 1) + (currentPE + 1 - remainderOfPermutation) * numberPermutationPEs;
		}
	}

	// int iter = 0;

	// compute only when some work is assigned to the current PE
	if (initial_perm_index <= end_perm_index) {
		vector<int> current_reply;

		vector<int>nodes_init = nth_permutation(nodes, initial_perm_index);
		vector<int>nodes_end = nth_permutation(nodes, end_perm_index);
		do
		{
			//Assign nodes to temporary
			vector<int>temp_Nodes = nodes_init;

			//Push last nodes is 0
			temp_Nodes.push_back(0);

			//Assign first nodes is 0
			temp_Nodes.insert(temp_Nodes.begin(), 0);

			/*for (int i = 0; i < temp.size(); i++) {
			printf("%d\n", temp[i]);
			}
			printf("\n");
			display_matrix(matrix);*/

			//currentPathCost is to get the current total path cost of the random matrix
			int currentPathCost = find_path_cost(matrix, temp_Nodes);

			/*printf("\n%d\n", val);*/

			//Check the lowest path cost then assign the cost path (currentPathCost)
			if (currentPathCost < optimal_cost)
			{
				optimal_cost = currentPathCost;
				current_reply = temp_Nodes;
			}

			// iter++;
			if (nodes_init == nodes_end) break;

		} while (next_permutation(nodes_init.begin(), nodes_init.end()));
		copy(current_reply.begin(), current_reply.end(), reply);

	}

	// comparing the optimal values from other PEs in the MASTERPE
	if (currentPE == MASTERPE) {

		// receiving reply from all other PEs
		for (int pe = 1; pe < totalPEs; pe++) {

			int tmp_optimal_value;
			int* tmp_ans = new int[N + 1];
			MPI_Recv(&tmp_ans[0], (N + 1), MPI_INT, pe, pe, MPI_COMM_WORLD, &stat);
			MPI_Recv(&tmp_optimal_value, 1, MPI_INT, pe, pe, MPI_COMM_WORLD, &stat);

			if (tmp_optimal_value < optimal_cost) {
				optimal_cost = tmp_optimal_value;
				copy(tmp_ans, tmp_ans + N + 1, reply);
			}
		}
	}
	else {
		// send reply to MASTERPE
		MPI_Send(&reply[0], (N + 1), MPI_INT, MASTERPE, currentPE, MPI_COMM_WORLD);
		MPI_Send(&optimal_cost, 1, MPI_INT, MASTERPE, currentPE, MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	endTime = MPI_Wtime();								// endTime time

	// Display the minimum cost path
	if (currentPE == MASTERPE) {
		printf("Minimum Cost Path for %d cities: \n", N);
		for (int i = 0; i < N + 1; i++) {
			cout << reply[i] << ' ';
		}
		cout << endl << endl;
	}

	// Display the minimum path cost
	if (currentPE == MASTERPE) {
		printf("Minimum Total Path Cost for %d cities: \n", N);
		int cost = 0;
		for (int i = 1; i < N + 1; i++)
		{
			cost += matrix[reply[i]][reply[i - 1]];
		}
		cout << cost << endl << endl;
	}

	// Display the run-time
	if (currentPE == MASTERPE) {
		printf("Run-Time for %d cities: ", N);
		cout << fixed << setprecision(5) << (endTime - startTime) << endl;
	}

	// deleting allocated memory
	delete reply;

	MPI_Finalize();

}