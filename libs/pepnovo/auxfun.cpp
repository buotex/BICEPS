#include "auxfun.h"

//#include "boost/filesystem/operations.hpp"
//#include "boost/filesystem/path.hpp"


int read_files_in_dir(const string& dir_path, vector<string>& file_names)
{

	return 0;
}



static unsigned int SEED;	

void rand_seed (unsigned int init)   {
	if (init != 0)
	{
		SEED = init;
	}
	else
	{
		time_t ltime;
		unsigned int t=(unsigned int)time( &ltime );

		SEED = t;
	}
}

unsigned int get_random_seed() { return SEED; }


/* Returns random uniform number */
double my_random()  
{
  static unsigned int a = 1588635695, m = 4294967291U, q = 2, r = 1117695901;

   SEED = a*(SEED % q) - r*(SEED / q);
   return ((double)SEED / (double)m);
}

void error()
{
	cout << "Error!!!" << endl;
	exit(1);
}

void error(const char *msg)
{
	cout << "Error: " << msg << endl;
	exit(1);
}

// returns the file size in bytes
int getFileSize(const char* sFileName)
{
  ifstream f;
  f.open(sFileName, ios_base::binary | ios_base::in);
  if (!f.good() || f.eof() || !f.is_open()) { return 0; }
  f.seekg(0, ios_base::beg);
  ifstream::pos_type begin_pos = f.tellg();
  f.seekg(0, ios_base::end);
  return static_cast<int>(f.tellg() - begin_pos);
}



void get_file_name_without_extension(const string& full_path, string& file_name)
{
	const int last_pos = full_path.length()-1;
	int dot_pos = last_pos;
	while (dot_pos>=0 && full_path[dot_pos] != '.')
		dot_pos--;

	if (dot_pos <=0)
		dot_pos = last_pos;

	int slash_pos=last_pos;
	while (slash_pos>=0 && (full_path[slash_pos] != '/') && (full_path[slash_pos] != '\\'))
		slash_pos--;

	file_name = full_path.substr(slash_pos+1,dot_pos-slash_pos-1);
}


// chooses k numbers from 0,...,n-1 (unique)
struct ChoosePair {
	ChoosePair(int i,double d) : idx(i),val(d) {};
	bool operator < (const ChoosePair& other) const
	{
		return val<other.val;
	}
	int idx;
	double val;
};

void choose_k_from_n(int k, int n, vector<int>& idxs)
{
	int i;
	idxs.clear();
	idxs.resize(k);

	vector<ChoosePair> pairs;

	if (k>n)
	{
		cout << "Error: choose " << k << " from " << n << " !" << endl;
		exit(1);
	}

	for (i=0; i<n; i++)
		pairs.push_back(ChoosePair(i,my_random()));
	
	sort(pairs.begin(),pairs.end());

	for (i=0; i<k; i++)
		idxs[i]=pairs[i].idx;

	sort(idxs.begin(),idxs.end());
}


mass_t ppm_val(mass_t offset, mass_t total_mass)
{
	return (offset / total_mass) * 1000000;
}



/*************************************************************
   finds all the permutaitons of n elements, repeated elements
   are allowed and do not create redundant permutations.
**************************************************************/
void generate_all_permutations(const vector<int>& org_vector, 
							   vector< vector<int> >& permutations)
{
	int i;
	vector<int> counts, symbols;
	permutations.clear();

	if (org_vector.size() == 0)
		return;

	counts.clear();
	symbols.clear();

	// create vector with symbols and their counts
	symbols.push_back(org_vector[0]);
	counts.push_back(1);

	for (i=1; i<org_vector.size(); i++)
	{
		int j;
		for (j=0; j<counts.size(); j++)
		{
			if (org_vector[i] == symbols[j])
			{
				counts[j]++;
				break;
			}
		}

		if (j == counts.size())
		{
			symbols.push_back(org_vector[i]);
			counts.push_back(1);
		}
	}

	vector<int> next_sym_idx,perm;
	int n = org_vector.size(); // total number of elements
	int k = counts.size(); // total number of element types
	next_sym_idx.resize(n,0);
	perm.resize(n,-1);
	int d=0;

	while (1)
	{
		while (next_sym_idx[d]<k && counts[next_sym_idx[d]] == 0)
			next_sym_idx[d]++;

		if (next_sym_idx[0]==k)
			break;

		if (next_sym_idx[d] >= k)
		{
			next_sym_idx[d]=0;
			d--;
			counts[next_sym_idx[d]]++;
			next_sym_idx[d]++;
			continue;
		}

		// add symbol
		perm[d]=symbols[next_sym_idx[d]];
		counts[next_sym_idx[d]]--;
		d++;

		if (d == n)
		{
			permutations.push_back(perm);
	//		int k;
	//		for (k=0; k<perm.size(); k++)
	//			cout << perm[k] << " ";
	//		cout << endl;

			d--;
			counts[next_sym_idx[d]]++;
			next_sym_idx[d]++;
		}
	}
}

// returns the minimal x for which the cumulative probability
// P(X<x)>= target_prob, assuming X~bin(n,p)
int get_min_number_from_binomial_prob(int n, double p, double target_prob)
{
	const double one_minus_p = 1.0 - p;

	double pow_p=1.0;
	double pow_1_minus_p = pow(one_minus_p,n);

	double sum_prob = pow_1_minus_p;
	double bin_coef = 1.0;
	double pow_val  = pow_1_minus_p;

//	cout << 0 << " " <<  pow_val << " " << pow_val << endl;
	int b=0;
	while (sum_prob<target_prob)
	{
		b++;
		bin_coef *= (double)(n-b+1);
		bin_coef /= (double)b;

		pow_val *= p;
		pow_val /= one_minus_p;

		double prob = bin_coef * pow_val;
		sum_prob += prob;

	//	cout << b << " " << prob << " " << sum_prob << endl;
	}
	return b;
}


void add_to_mass_vector(vector<mass_t>& vec, mass_t val, mass_t tolerance)
{
	int i;
	for (i=0; i<vec.size(); i++)
		if (fabs(vec[i]-val)<tolerance)
			break;
	if (i<vec.size())
		return;
	vec.push_back(val);
	
}


/******************************************************************************
splits string according to a given delimeter char.
*******************************************************************************/
void split_string(const string& str, vector<string>& results, char delim)
{
	const char *str_buff = str.c_str();
	const int   max_pos  = str.length();

	results.clear();
	int last_pos = 0;
	int pos=0;

	for (pos=0; pos<max_pos; pos++)
	{
		if (str_buff[pos]==delim) 
		{
			if (pos>last_pos)
				results.push_back(str.substr(last_pos,pos-last_pos));
			last_pos=pos+1;
			continue;
		}
		else if (pos == max_pos-1)
		{
			if (pos>last_pos)
				results.push_back(str.substr(last_pos,pos-last_pos+1));
			break;
		}
	}
}
