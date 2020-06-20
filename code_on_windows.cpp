#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <utility>
#include <math.h>
#include <algorithm>
#include <tuple>
using namespace std;

using doublepair = pair<double, double>;
using cliquetype = std::tuple<vector<doublepair>, vector<vector<int>>>;

bool removable(vector<int> &neighbor, vector<int> &cover);
int max_removable(vector<vector<int> > &neighbors, vector<int> &cover);
vector<int> procedure_1(vector<vector<int>> &neighbors, vector<int> &cover);
vector<int> procedure_2(vector<vector<int>> &neighbors, vector<int> &cover, int k);
int cover_size(vector<int> &cover);
vector<int> minCliqueCover(vector<vector<int>> &graph, int n);
vector<doublepair> minCliqueCover(vector<doublepair> &);
bool emptyGraph(vector<vector<int>> &graph);
void removeNode(vector<vector<int>> &graph, vector<int> &nodes);

bool emptyGraph(vector<vector<int>> &graph)
{
	bool empty = true;
	for (int i = 0; i < graph.size(); ++i)
	{
		if (std::find(graph[i].begin(), graph[i].end(), 0) != graph[i].end())
		{
			empty = false;
		}
	}
	return empty;
}

void removeNode(vector<vector<int>> &graph, vector<int> &nodes)
{
	for (int k = 0; k < nodes.size(); ++k)
	{
		for (int i = 0; i < graph.size(); ++i)
		{
			for (int j = 0; j < graph.size(); ++j)
			{
				if (i == nodes[k] || j == nodes[k])
				{
					graph[i][j] = 1;
				}
			}
		}
	}
}

bool intersect(const doublepair &set1, const doublepair &set2)
{
	double epsilon = 1e-3,
		input_min = set1.first, input_max = set1.second,
		output_min = set2.first - epsilon, output_max = set2.second + epsilon,
		step = 0.01;
	// set1 is input, set2 is output
	bool indicator = false;
	for (auto val = input_min; val <= input_max; val = val + step)
	{
		if (val >= output_min && val <= output_max)
		{

			indicator = true;
			break;
		}
	}
	return indicator;
}

vector<doublepair> minCliqueCover(vector<doublepair> &ranges)
{
	vector<int> all_index;
	cliquetype results;
	vector<doublepair> output_ranges;
	vector<vector<int>> graph(ranges.size()), indices;
	auto gbeg = graph.begin();
	for (int i = 0; i < ranges.size(); ++i)
	{
		vector<int> row(ranges.size());
		auto rbeg = row.begin();
		for (int j = 0; j < ranges.size(); ++j)
		{
			if (i == j)
			{
				*rbeg++ = 1;
			}
			else
			{
				if (intersect(ranges[i], ranges[j]))
				{
					*rbeg++ = 0;
				}
				else
				{
					*rbeg++ = 1;
				}
			}
		}
		*gbeg++ = row;
	}

	vector<int> index;
	while (!emptyGraph(graph))
	{
		for (int i = 0; i < ranges.size(); ++i)
		{
			for (int j = 0; j < ranges.size(); ++j)
			{
				cout << graph[i][j] << " ";
			}
			cout << endl;
		}

		cout << "graph built" << endl;

		index = minCliqueCover(graph, graph.size() - index.size());
		all_index.insert(all_index.end(), index.begin(), index.end());
		indices.push_back(index);
		for (int i = 0; i < index.size(); ++i)
		{
			cout << index[i] << ", " << endl;
		}
		vector<double> minvec, maxvec;
		for (int i = 0; i < index.size(); ++i)
		{
			minvec.push_back(ranges[index[i]].first);
			maxvec.push_back(ranges[index[i]].second);
		}
		output_ranges.push_back(std::make_pair(*std::max_element(minvec.begin(), minvec.end()),
			*std::min_element(maxvec.begin(), maxvec.end())));

		removeNode(graph, index);
	}
	for (int i = 0; i < ranges.size(); ++i)
	{
		if (std::find(all_index.begin(), all_index.end(), i) == all_index.end())
		{
			indices.insert(indices.end(), vector<int>{i});
			output_ranges.insert(output_ranges.end(), ranges[i]);
		}
	}
	results = std::make_tuple(output_ranges, indices);
	return output_ranges;
}


vector<int> minCliqueCover(vector<vector<int>> &graph, int _n)
{
	int n = graph.size(), K, i, j, k, p, q, r, s, min, edge;

	vector<vector<int>> neighbors;
	for (int i = 0; i < graph.size(); i++)
	{
		vector<int> neighbor;
		for (int j = 0; j < graph[i].size(); j++)
			if (graph[i][j] == 1)
				neighbor.push_back(j);
		neighbors.push_back(neighbor);
	}

	vector<int> clique_index;
	for (size_t K = _n; K != 1; K--)
	{
		std::cout << "K: " << K << std::endl;
		int counter = 0;

		k = n - K;
		//Find Cliques
		bool found = false;
		
		min = n + 1;
		vector<vector<int>> covers;
		vector<int> allcover(graph.size());
		for (i = 0; i < graph.size(); i++)
			allcover[i] = 1;
		for (i = 0; i < allcover.size(); i++)
		{
			if (found)
				break;
			counter++;
			vector<int> cover = allcover;
			cover[i] = 0;
			cover = procedure_1(neighbors, cover);
			s = cover_size(cover);
			if (s < min)
				min = s;
			if (s <= k)
			{
				if (n - s == K)
				{
					clique_index.resize(K);
					auto cibeg = clique_index.begin();
					for (j = 0; j < cover.size(); j++)
						if (cover[j] == 0)
							*cibeg++ = j;
				}
				covers.insert(covers.end(), cover);
				found = true;
				break;
			}
			for (j = 0; j < n - k; j++)
				cover = procedure_2(neighbors, cover, j);
			s = cover_size(cover);
			if (s < min)
				min = s;
			if (n-s == K)
			{
				clique_index.resize(K);
				auto cibeg = clique_index.begin();
				for (j = 0; j < cover.size(); j++)
					if (cover[j] == 0)
						*cibeg++ = j;
			}
			covers.insert(covers.end(), cover);
			if (s <= k)
			{
				found = true;
				break;
			}
		}
		//Pairwise Intersections
		for (p = 0; p < covers.size(); p++)
		{
			if (found)
				break;
			for (q = p + 1; q < covers.size(); q++)
			{
				if (found)
					break;
				counter++;
				
				vector<int> cover = allcover;
				for (r = 0; r < cover.size(); r++)
					if (covers[p][r] == 0 && covers[q][r] == 0)
						cover[r] = 0;
				cover = procedure_1(neighbors, cover);
				s = cover_size(cover);
				if (s < min)
					min = s;
				if (s <= k)
				{
					if (n-s ==K)
					{
						clique_index.resize(K);
						auto cibeg = clique_index.begin();
						for (j = 0; j < cover.size(); j++)
							if (cover[j] == 0)
								*cibeg++ = j;
					}
					found = true;
					break;
				}
				for (j = 0; j < k; j++)
					cover = procedure_2(neighbors, cover, j);
				s = cover_size(cover);
				if (s < min)
					min = s;
				if (n - s == K)
				{
					clique_index.resize(K);
					auto cibeg = clique_index.begin();
					for (j = 0; j < cover.size(); j++)
						if (cover[j] == 0)
							*cibeg++ = j;
				}
				if (s <= k)
				{
					found = true;
					break;
				}
			}
		}
		if (found)
			break;
	}
	return clique_index;	
}



int main()
{

	vector<pair<double, double>> ranges{
		{6.0, 15.0}, {4.5, 4.5}, {6.0, 15.0}, {4.5, 4.5}, {6.0, 15.0}, {4.5, 4.5}, {6.0, 15.0}, {4.5, 4.5},
		{6.0, 15.0}, {4.5, 4.5}, {6.0, 15.0}, {4.5, 4.5}, {6.0, 15.0}, {4.5, 4.5}, {6.0, 15.0}, {4.5, 4.5},
		{6.0, 15.0}, {4.5, 4.5}, {6.0, 15.0}, {4.5, 4.5}, {6.0, 15.0}, {4.5, 4.5}, {6.0, 15.0}, {4.5, 4.5},
		{6.0, 15.0}, {4.5, 4.5}, {6.0, 15.0}, {4.5, 4.5}, {6.0, 15.0}, {4.5, 4.5}, {6.0, 15.0}, {4.5, 4.5},
//		{4.5, 4.5} , {6.0, 9.0}, {6.0, 8.0}, {12.0, 15.0}
	};
	
	vector<doublepair> output_ranges = minCliqueCover(ranges);
	for (int i = 0; i < output_ranges.size(); ++i)
	{
		cout << output_ranges[i].first << ", " << output_ranges[i].second << std::endl;
	}

/*
	vector<vector<int>> graph(ranges.size());
	auto gbeg = graph.begin();
	for (int i = 0; i < ranges.size(); ++i)
	{
		vector<int> row(ranges.size());
		auto rbeg = row.begin();
		for (int j = 0; j < ranges.size(); ++j)
		{
			if (i == j)
			{
				*rbeg++ = 1;
			}
			else
			{
				if (intersect(ranges[i], ranges[j]))
				{
					*rbeg++ = 0;
				}
				else
				{
					*rbeg++ = 1;
				}
			}
		}
		*gbeg++ = row;
	}
	for (int i = 0; i < ranges.size(); ++i)
	{
		for (int j = 0; j < ranges.size(); ++j)
		{
			cout << graph[i][j] << " ";
		}
		cout << endl;
	}
	cout << "graph built" << endl;
	auto index = minCliqueCover(graph);
	for (int i = 0; i < index.size(); ++i)
	{
		cout << index[i] << ", " << ends;
	}
	system("pause");
	return 0;
	
	/*
	//Read Graph (note we work with the complement of the input graph)
	cout << "Clique Algorithm." << endl;
	int n, i, j, k, K, p, q, r, s, min, edge, counter = 0;
	infile >> n;
	vector<vector<int>> graph;
	for (i = 0; i < n; i++)
	{
		vector<int> row;
		for (j = 0; j < n; j++)
		{
			infile >> edge;
			if (edge == 0)
				row.push_back(1);
			else
				row.push_back(0);
		}
		graph.push_back(row);
	}
	//Find Neighbors
	vector<vector<int> > neighbors;
	for (i = 0; i < graph.size(); i++)
	{
		vector<int> neighbor;
		for (j = 0; j < graph[i].size(); j++)
			if (graph[i][j] == 1)
				neighbor.push_back(j);
		neighbors.push_back(neighbor);
	}
	cout << "Graph has n = " << n << " vertices." << endl;
	//Read maximum size of Clique wanted
	cout << "Find a Clique of size at least k = ";
	cin >> K;
	k = n - K;
	//Find Cliques
	bool found = false;
	cout << "Finding Cliques..." << endl;
	min = n + 1;
	vector<vector<int> > covers;
	vector<int> allcover;
	for (i = 0; i < graph.size(); i++)
		allcover.push_back(1);
	for (i = 0; i < allcover.size(); i++)
	{
		if (found)
			break;
		counter++;
		cout << counter << ". ";
		outfile << counter << ". ";
		vector<int> cover = allcover;
		cover[i] = 0;
		cover = procedure_1(neighbors, cover);
		s = cover_size(cover);
		if (s < min)
			min = s;
		if (s <= k)
		{
			outfile << "Clique (" << n - s << "): ";
			for (j = 0; j < cover.size(); j++)
				if (cover[j] == 0)
					outfile << j + 1 << " ";
			outfile << endl;
			cout << "Clique Size: " << n - s << endl;
			covers.push_back(cover);
			found = true;
			break;
		}
		for (j = 0; j < n - k; j++)
			cover = procedure_2(neighbors, cover, j);
		s = cover_size(cover);
		if (s < min)
			min = s;
		outfile << "Clique (" << n - s << "): ";
		for (j = 0; j < cover.size(); j++)
			if (cover[j] == 0)
				outfile << j + 1 << " ";
		outfile << endl;
		cout << "Clique Size: " << n - s << endl;
		covers.push_back(cover);
		if (s <= k)
		{
			found = true;
			break;
		}
	}
	//Pairwise Intersections
	for (p = 0; p < covers.size(); p++)
	{
		if (found)
			break;
		for (q = p + 1; q < covers.size(); q++)
		{
			if (found)
				break;
			counter++;
			cout << counter << ". ";
			outfile << counter << ". ";
			vector<int> cover = allcover;
			for (r = 0; r < cover.size(); r++)
				if (covers[p][r] == 0 && covers[q][r] == 0)
					cover[r] = 0;
			cover = procedure_1(neighbors, cover);
			s = cover_size(cover);
			if (s < min)
				min = s;
			if (s <= k)
			{
				outfile << "Clique (" << n - s << "): ";
				for (j = 0; j < cover.size(); j++)
					if (cover[j] == 0)
						outfile << j + 1 << " ";
				outfile << endl;
				cout << "Clique Size: " << n - s << endl;
				found = true;
				break;
			}
			for (j = 0; j < k; j++)
				cover = procedure_2(neighbors, cover, j);
			s = cover_size(cover);
			if (s < min)
				min = s;
			outfile << "Clique (" << n - s << "): ";
			for (j = 0; j < cover.size(); j++)
				if (cover[j] == 0)
					outfile << j + 1 << " ";
			outfile << endl;
			cout << "Clique Size: " << n - s << endl;
			if (s <= k)
			{
				found = true;
				break;
			}
		}
	}
	if (found)
		cout << "Found Clique of size at least " << K << "." << endl;
	else
		cout << "Could not find Clique of size at least " << K << "." << endl
				<< "Maximum Clique size found is " << n - min << "." << endl;
	cout << "See cliques.txt for results." << endl;
	return 0;
	*/

}

bool removable(vector<int> &neighbor, vector<int> &cover)
{
	bool check = true;
	for (int i = 0; i < neighbor.size(); i++)
		if (cover[neighbor[i]] == 0)
		{
			check = false;
			break;
		}
	return check;
}

int max_removable(vector<vector<int>> &neighbors, vector<int> &cover)
{
	int r = -1, max = -1;
	for (int i = 0; i < cover.size(); i++)
	{
		if (cover[i] == 1 && removable(neighbors[i], cover) == true)
		{
			vector<int> temp_cover = cover;
			temp_cover[i] = 0;
			int sum = 0;
			for (int j = 0; j < temp_cover.size(); j++)
				if (temp_cover[j] == 1 && removable(neighbors[j], temp_cover)
					== true)
					sum++;
			if (sum > max)
			{
				max = sum;
				r = i;
			}
		}
	}
	return r;
}

vector<int> procedure_1(vector<vector<int> > &neighbors, vector<int> &cover)
{
	vector<int> temp_cover = cover;
	int r = 0;
	while (r != -1)
	{
		r = max_removable(neighbors, temp_cover);
		if (r != -1)
			temp_cover[r] = 0;
	}
	return temp_cover;
}

vector<int> procedure_2(vector<vector<int> > &neighbors, vector<int> &cover,
	int k)
{
	int count = 0;
	vector<int> temp_cover = cover;
	int i = 0;
	for (int i = 0; i < temp_cover.size(); i++)
	{
		if (temp_cover[i] == 1)
		{
			int sum = 0, index;
			for (int j = 0; j < neighbors[i].size(); j++)
				if (temp_cover[neighbors[i][j]] == 0)
				{
					index = j;
					sum++;
				}
			if (sum == 1 && cover[neighbors[i][index]] == 0)
			{
				temp_cover[neighbors[i][index]] = 1;
				temp_cover[i] = 0;
				temp_cover = procedure_1(neighbors, temp_cover);
				count++;
			}
			if (count > k)
				break;
		}
	}
	return temp_cover;
}

int cover_size(vector<int>& cover)
{
	int count = 0;
	for (int i = 0; i < cover.size(); i++)
		if (cover[i] == 1)
			count++;
	return count;
}

