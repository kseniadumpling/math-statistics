#include <algorithm>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
#include <iomanip>

using namespace std;

vector<double> normal(int count, double mean = 0, double std_dev = 1) {
	vector<double> vec(count);
	random_device randomness_device{};
	mt19937 pseudorandom_generator{ randomness_device() };
	normal_distribution<double> distribution{ mean, std_dev };
	for (int i = 0; i < count; i++)
	{
		vec[i] = distribution(pseudorandom_generator);
	}

	return vec;
}

vector<double> mixture(int count, vector<double> mean = { 0 }, vector<double> std_dev = { 1 }, initializer_list<double> weight = { 1 }) {
	vector<double> vec(count);
	random_device randomness_device{};
	mt19937 pseudorandom_generator{ randomness_device() };
	vector<normal_distribution<double>> distribution_normal(mean.size());
	for (int i = 0; i< distribution_normal.size() ;i++)
	{
		distribution_normal[i] = normal_distribution<double>{ mean[i], std_dev[i] };
	}
	discrete_distribution<int> w{ weight };
	int index;
	for (int i = 0; i < count; i++)
	{
		index = w(pseudorandom_generator);
		vec[i] = distribution_normal[index](pseudorandom_generator);
	}

	return vec;
}

vector<double> cauchy(int count, double mean = 0, double std_dev = 1) {
	vector<double> vec(count);
	random_device randomness_device{};
	mt19937 pseudorandom_generator{ randomness_device() };
	cauchy_distribution<double> distribution{ mean, std_dev };
	for (int i = 0; i < count; i++)
	{
		vec[i] = distribution(pseudorandom_generator);
	}

	return vec;
}

vector<double> laplas(int count, double mean = 0, double std_dev = 1) {
	vector<double> vec(count);
	random_device randomness_device{};
	mt19937 pseudorandom_generator{ randomness_device() };
	exponential_distribution<double> distribution{ 1.0 / std_dev };
	discrete_distribution<int> w{ {0.5,0.5} };
	for (int i = 0; i < count; i++)
	{
		vec[i] = (w(pseudorandom_generator) == 1 ? 1 : -1) * distribution(pseudorandom_generator) + mean;
	}

	return vec;
}

vector<double> uniform(int count, double min = 0, double max = 1) {
	vector<double> vec(count);
	random_device randomness_device{};
	mt19937 pseudorandom_generator{ randomness_device() };
	uniform_real<double> distribution{ min, max };
	for (int i = 0; i < count; i++)
	{
		vec[i] = distribution(pseudorandom_generator);
	}

	return vec;
}

double sampleMean(vector<double> vec) {
	int n = vec.size();
	double a = 0.0;
	for (double item : vec)
	{
		a += item;
	}
	return a / n;
}

double sampleMeanSquare(vector<double> vec) {
	int n = vec.size();
	double a = 0.0;
	for (double item : vec)
	{
		a += item * item;
	}
	return a / n;
}

double Dz(vector<double> vec) {
	double sm = sampleMean(vec);
	return sampleMeanSquare(vec) - sm * sm;
}

double median(vector<double> vec) {
	int n = vec.size();
	if (n % 2 == 0) {
		return (vec[n / 2 - 1] + vec[n / 2]) / 2;
	}
	else {
		return vec[(int)(n / 2)];
	}
}

double extrimeMean(vector<double> vec) {
	int n = vec.size();
	return (vec[n-1] + vec[0]) / 2;
}

double trimmedMean(vector<double> vec, double r = 0.05) {
	int n = vec.size();
	double a = 0.0;
	for (int i = r * n; i < n * (1.0 - r); i++)
	{
		a += vec[i];
	}
	return a / (n - (int) 2 * r * n);
}

void print_result(ofstream &fout, const char* str, vector<double> a, vector<double> b, vector<double> c, vector<double> d) {
	fout << str << endl <<
		'\t' << setw(15) << 'a' << setw(15) << 'b' << setw(15) << 'c' << setw(15) << 'd' << endl
		<< "Zm" << '\t'
		<< setw(15) << sampleMean(a)
		<< setw(15) << sampleMean(b)
		<< setw(15) << sampleMean(c)
		<< setw(15) << sampleMean(d) << endl
		<< "Zm2" << '\t'
		<< setw(15) << sampleMeanSquare(a)
		<< setw(15) << sampleMeanSquare(b)
		<< setw(15) << sampleMeanSquare(c)
		<< setw(15) << sampleMeanSquare(d) << endl
		<< "Dz" << '\t'
		<< setw(15) << Dz(a)
		<< setw(15) << Dz(b)
		<< setw(15) << Dz(c)
		<< setw(15) << Dz(d) << endl << endl;

	return;
}

int main() {
	int n = 100;
	double r = 0.05;
	int m = 1000;
	// normal(n);
	// uniform(n, -sqrt(3), sqrt(3));
	// laplas(n, 0, sqrt(2) / 2);
	// cauchy(n);
	// mixture(n, { 0.0, 0.0 }, { 1.0, 3.0 }, { 0.9, 0.1 });
	ofstream fout;
	fout.open("test.txt");
	vector<double> vec(n);
	vector<double> a(m);
	vector<double> b(m);
	vector<double> c(m);
	vector<double> d(m);

	for (int i = 0; i < m; i++) {
		vec = normal(n);
		sort(vec.begin(), vec.end());
		a[i] = sampleMean(vec);
		b[i] = median(vec);
		c[i] = extrimeMean(vec);
		d[i] = trimmedMean(vec, r);
	}
	print_result(fout, "normal", a, b, c, d);

	for (int i = 0; i < m; i++) {
		vec = uniform(n, -sqrt(3), sqrt(3));
		sort(vec.begin(), vec.end());
		a[i] = sampleMean(vec);
		b[i] = median(vec);
		c[i] = extrimeMean(vec);
		d[i] = trimmedMean(vec, r);
	}
	print_result(fout, "uniform", a, b, c, d);

	for (int i = 0; i < m; i++) {
		vec = laplas(n, 0, sqrt(2) / 2);
		sort(vec.begin(), vec.end());
		a[i] = sampleMean(vec);
		b[i] = median(vec);
		c[i] = extrimeMean(vec);
		d[i] = trimmedMean(vec, r);
	}
	print_result(fout, "laplas", a, b, c, d);

	for (int i = 0; i < m; i++) {
		vec = cauchy(n);
		sort(vec.begin(), vec.end());
		a[i] = sampleMean(vec);
		b[i] = median(vec);
		c[i] = extrimeMean(vec);
		d[i] = trimmedMean(vec, r);
	}
	print_result(fout, "cauchy", a, b, c, d);

	for (int i = 0; i < m; i++) {
		vec = mixture(n, { 0.0, 0.0 }, { 1.0, 3.0 }, { 0.9, 0.1 });
		sort(vec.begin(), vec.end());
		a[i] = sampleMean(vec);
		b[i] = median(vec);
		c[i] = extrimeMean(vec);
		d[i] = trimmedMean(vec, r);
	}
	print_result(fout, "mixture", a, b, c, d);

	fout.close();

	return 0;
}