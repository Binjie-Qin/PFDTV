// VS2015+opencv3.1+matlab2013
#include <iostream>
#include "engine.h"
#include<opencv2\opencv.hpp>
#include<opencv2\highgui.hpp>
#include<stdio.h>
#include<cmath>
#include<fstream>
#include<time.h>


#pragma comment(lib, "libeng.lib")
#pragma comment(lib, "libmx.lib")
#pragma comment(lib, "libmex.lib")

using namespace std;
using namespace cv;



Mat calorder(Mat inputarray)
{
	int M = inputarray.rows;
	int N = inputarray.cols;
	float x = 0;
	Mat_<float> ordernor(M, N);
	Mat_<float> order(M, N);

	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			ordernor.at<float>(i, j) = inputarray.at<float>(i, j)*inputarray.at<float>(i, j);
			//x = 1 + ordernor.at<float>(i, j) / 1.44f;
			x = 1 + ordernor.at<float>(i, j);
			order.at<float>(i, j) = 1 + log(x) / log(2);
		}
	}
	return order;
}

Mat normal(Mat input)
{
	Mat_<float> output(input.rows, input.cols);
	float min = input.at<float>(0, 0);
	float max = input.at<float>(0, 0);
	for (int i = 0; i < input.rows; i++)
		for (int j = 0; j < input.cols; j++)
		{
			if (input.at<float>(i, j) > max)
				max = input.at<float>(i, j);
			if (input.at<float>(i, j) < min)
				min = input.at<float>(i, j);
		}
	for (int i = 0; i < input.rows; i++)
		for (int j = 0; j < input.cols; j++)
			output.at<float>(i, j) = (input.at<float>(i, j) - min) / (max - min);
	return output;
}

// the first weighted coefficient
Mat xishu1(Mat input)
{
	int M = input.rows;
	int N = input.cols;
	Mat_<float>output(M, N);
	
	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
			output.at<float>(i, j) = input.at<float>(i, j)*(2 - input.at<float>(i, j));
	return output;
}

// the second weighted coefficient
Mat xishu2(Mat input)
{
	int M = input.rows;
	int N = input.cols;
	Mat_<float>output(M, N);
	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
			output.at<float>(i, j) = pow(input.at<float>(i, j) - 1, 2);
	return output;
}

Mat AFGOM(engine *ep, Mat gray, Mat gray0)
{
	int k0 = 20;
	float k1;
	k1 = k0*exp(-0.05*(8-1));
	int M = gray.rows;
	int N = gray.cols;
	Mat_<float>TVDX(M, N);//FTVX
	Mat_<float>TVDY(M, N);//FTVY
	Mat_<float>PMDX(M, N);//FADX 
	Mat_<float>PMDY(M, N);//FADY
	Mat_<float>order(M, N);
	Mat_<float>xi1(M, N);
	Mat_<float>xi2(M, N);
	Mat_<float>U(M, N);
	Mat_<float>PCMAP(M, N);
	
	// calculate the phase asymmetry by using matlab script
	// matlab script path
	engEvalString(ep, "addpath('C:/Users/HD/Documents/Visual Studio 2015/Projects/ConsoleApplication6/ConsoleApplication6');");
	//engEvalString(ep, "test25;");
	engEvalString(ep, "PAS;");


	// the phase asymmetry result is saved in a txt
	// generate the new txt after each iteration
	fstream file;
	file.open("D:\\matlab\\PCMAP.txt");
	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
			file >> PCMAP.at<float>(i, j);
	file.close();


	PCMAP = normal(PCMAP);
	order = calorder(PCMAP);

	Mat_<float> si(M, N);
	Mat_<float> ita(M, N);
	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
		{
			//PCMAP.at<float>(i, j) = pow(PCMAP.at<float>(i, j), 2);
			si.at<float>(i, j) = 1 + PCMAP.at<float>(i, j) * 254;
		}
	
	xi1 = xishu1(PCMAP);
	xi2 = xishu2(PCMAP);
	  
	float test;
	float result;
	float middle;

	Mat_<float>delta_x(M, N);
	Mat_<float>delta_y(M, N);

	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
			for (int k = 0; k < j + 1; k++)
			{
				test = tgammaf(order.at<float>(i, j) + 1) / (tgammaf(k + 1)*tgammaf(order.at<float>(i, j ) - k + 1));
				if (cvIsNaN(test) == 1)
					test = 0;
				result = delta_x.at<float>(i, j) + pow(-1, k)*test*gray.at<float>(i, j - k);
				delta_x.at<float>(i, j) = result;
			}

	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
			for (int k = 0; k < i + 1; k++)
			{
				test = tgammaf(order.at<float>(i , j) + 1) / (tgammaf(k + 1)*tgammaf(order.at<float>(i , j) - k + 1));
				if (cvIsNaN(test) == 1)
					test = 0;
				result = delta_y.at<float>(i, j) + pow(-1, k)*test*gray.at<float>(i - k, j);
				delta_y.at<float>(i, j) = result;
			}

	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
			for (int l = 0; l < N - j; l++)
			{
				test = tgammaf(order.at<float>(i, j ) + 1) / (tgammaf(l + 1)*tgammaf(order.at<float>(i, j ) - l + 1));
				if (cvIsNaN(test) == 1)
					test = 0;
				result = TVDX.at<float>(i, j) + pow(-1, l)*test*delta_x.at<float>(i, j + l) / sqrt(pow(delta_x.at<float>(i, j + l), 2) + pow(delta_y.at<float>(i, j + l), 2) + 0.0001);
				TVDX.at<float>(i, j) = result;
				middle = PMDX.at<float>(i, j) + pow(-1, l)*test*pow(k1, 2)*delta_x.at<float>(i, j + l) / (k1*k1 + (pow(delta_x.at<float>(i, j + l), 2) + pow(delta_y.at<float>(i, j + l), 2))*pow(si.at<float>(i,j+l ),2)+ 0.0001);
				PMDX.at<float>(i, j) = middle;
			}
	
	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
			for (int l = 0; l < M - i; l++)
			{
				test = tgammaf(order.at<float>(i, j) + 1) / (tgammaf(l + 1)*tgammaf(order.at<float>(i, j) - l + 1));
				if (cvIsNaN(test) == 1)
					test = 0;
				result = TVDY.at<float>(i, j) + pow(-1, l)*test*delta_y.at<float>(i + l, j) / sqrt(pow(delta_x.at<float>(i + l, j), 2) + pow(delta_y.at<float>(i + l, j), 2) + 0.0001);
				TVDY.at<float>(i, j) = result;
				middle = PMDY.at<float>(i, j) + pow(-1, l)*test*pow(k1, 2)*delta_y.at<float>(i + l, j) / (k1*k1 + (pow(delta_x.at<float>(i + l, j), 2) + pow(delta_y.at<float>(i + l, j), 2))*pow(si.at<float>(i+l ,j),2) + 0.0001);
				PMDY.at<float>(i, j) = middle;
			}

	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
		{
			U.at<float>(i, j) = xi1.at<float>(i, j)*(TVDX.at<float>(i, j) + TVDY.at<float>(i, j)) + xi2.at<float>(i, j)*(PMDX.at<float>(i, j) + PMDY.at<float>(i, j));
			result = gray.at<float>(i, j) - 0.15*(U.at<float>(i, j) + 0.01*(gray.at<float>(i, j) - gray0.at<float>(i, j)));
			gray.at<float>(i, j) = result;
		}
	return gray;
}

int main()
{	
	clock_t start, finish;
	double totaltime;
	start = clock();  
	
	Engine* ep;
	ep = engOpen(NULL);

	engSetVisible(ep, false);
	if (ep == NULL)
	{
		cout << "Can't start MATLAB engine!";
		exit(-1);
	}
	string st = "./"; 

	Mat picture = imread("C:\\Users\\HD\\Desktop\\origin.jpg");
	Mat gray;
	cvtColor(picture, gray, CV_RGB2GRAY);
	gray.convertTo(gray, CV_32FC1);

	int M, N;
	M = gray.rows;
	N = gray.cols;

	// imread as u0£¬picture and picture0 have the same image path
	Mat picture0 = imread("C:\\Users\\HD\\Desktop\\origin.jpg");
	Mat gray0;
	cvtColor(picture0, gray0, CV_RGB2GRAY);
	gray0.convertTo(gray0, CV_32FC1);


	int n = 1;
	int Nit = 8;
	Mat_<float> tupian(M, N);

	while (n<Nit)
	{

		cout << "running"<< n << "iteration..." << endl;
		gray = AFGOM(ep, gray, gray0);
		

		for (int i = 0; i < M; i++)
			for (int j = 0; j < N; j++)
				tupian.at<float>(i, j) = gray.at<float>(i, j);
		gray.convertTo(gray, CV_8UC1);
		imwrite(st + "origin.jpg", gray); 
		gray.convertTo(gray, CV_32FC1);
		for (int i = 0; i < M; i++)
			for (int j = 0; j < N; j++)
				gray.at<float>(i, j) = tupian.at<float>(i, j);
		n++;
	}
	gray.convertTo(gray, CV_8UC1);

	imshow("test", gray);
	imwrite(st + "1.jpg", gray);
		
	engClose(ep);

	finish = clock();
	totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
	cout << "\n running time is: " << totaltime << "s" << endl;

	waitKey(2017041400);
}



