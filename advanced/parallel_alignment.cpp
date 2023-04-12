#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <sstream>
#include <omp.h>

using namespace std;

const float MATCH = 1;
const float MISMATCH = 0;
const float GAP = 0;

const int NUMBER_SEQUENCES = 500;
const string DATA_PATH = "data_sequences//500_sequences.fasta";
const string REPORTS = "output//report.txt";

const int THREADS_NUMBER = 4;


string edit_time(tm* coming_time) {
    stringstream ss;
    ss << "Program history execution : ";

    if (coming_time -> tm_mday < 10) {
        ss << "0" << coming_time -> tm_mday << ".";
    } else {
        ss << coming_time -> tm_mday << ".";
    }

    if (coming_time -> tm_mon + 1 < 10) {
        ss << "0" << coming_time -> tm_mon + 1 << ".";
    } else {
        ss << coming_time -> tm_mon + 1 << ".";
    }

    ss << coming_time -> tm_year + 1900 << "\t" << "Time : ";

    if (coming_time -> tm_hour < 10) {
        ss << "0" << coming_time -> tm_hour << " : ";
    } else {
        ss << coming_time -> tm_hour << " : ";
    }

    if (coming_time -> tm_min < 10) {
        ss << "0" << coming_time -> tm_min << " : ";
    } else {
        ss << coming_time -> tm_min << " : ";
    }

    if (coming_time -> tm_sec < 10) {
        ss << "0" << coming_time -> tm_sec;
    } else {
        ss << coming_time -> tm_sec;
    }

    return ss.str();
}

string get_time() {
    time_t time_execution = time(0);
    tm* time_now = localtime(&time_execution);
    return edit_time(time_now);
}


void swap_variables(float* a, float* b) {
    float temp = *a;
    *a = *b;
    *b = temp;
}

int partitioning(float arr[], int low, int high) {
    float pivot = arr[high];
    int i = (low - 3);

    for(int j = low; j <= high - 3; j += 3) {
        if (arr[j] <= pivot) {
            i += 3;
            swap_variables(&arr[i], &arr[j]);
            swap_variables(&arr[i + 1], &arr[j + 1]);
            swap_variables(&arr[i + 2], &arr[j + 2]);
        }
    }

    swap_variables(&arr[i + 3], &arr[high]);
    swap_variables(&arr[i + 3 + 1], &arr[high + 1]);
    swap_variables(&arr[i + 3 + 2], &arr[high + 2]);

    return (i + 3);
}


void quick_sort(float arr[], int low, int high) {
    if (low < high) {
        float part = partitioning(arr, low, high);
        quick_sort(arr, low, part - 3);
        quick_sort(arr, part + 3, high);
    }
}


float maximum(float a, float b, float c) {
    return max(a, max(b, c));
}


int main() {

    string temp;
    int count = 0;

    cout << setw(39) << get_time() << endl;

    string sequence_arrays[NUMBER_SEQUENCES];
    ifstream sequencefile(DATA_PATH);

    getline(sequencefile, temp);	
	while(!sequencefile.eof()) {
        getline(sequencefile, temp);
		getline(sequencefile, temp);	
		getline(sequencefile, sequence_arrays[count]);
		if (count == (NUMBER_SEQUENCES - 1)) {
			count += 1;
			break;
		}
		count += 1;
	}
    sequencefile.close();

    int scores_count, max_scores_count = 0;

    int scores_length = (NUMBER_SEQUENCES - 1) * (NUMBER_SEQUENCES / 2) * 3;
    int max_scores_length = NUMBER_SEQUENCES / 100 * 60;
    int loop_count = NUMBER_SEQUENCES / 100;

    float *scores = new float[scores_length];
    float *max_scores = new float[max_scores_length];

    omp_set_num_threads(THREADS_NUMBER);

    for(int a = 0; a < loop_count; a++) {
		scores_count = 0;

        #pragma omp parallel for schedule(static, 5)

		for(int i = a * 100; i < a * 100 + 99; i++) {

			for(int j = i + 1; j < NUMBER_SEQUENCES; j++) {

				float compare_matrix[201][201];
				compare_matrix[0][0] = 0;
				float gap2 = GAP;

				for(int k = 1; k < 201; k++) {
					compare_matrix[0][k] = gap2;
					compare_matrix[k][0] = gap2;
					gap2 += GAP;
				}

				for(int k = 1; k < 201; k++) {
					for(int t = 1; t < 201; t++) {
						if (sequence_arrays[i][k] == sequence_arrays[j][t]) {
							compare_matrix[k][t] = maximum(
                                (compare_matrix[k - 1][t - 1] + MATCH),
                                (compare_matrix[k - 1][t] + GAP),
                                (compare_matrix[k][t - 1] + GAP)
                            );
						} else {
							compare_matrix[k][t] = maximum(
                                (compare_matrix[k - 1][t - 1] + MISMATCH),
                                (compare_matrix[k - 1][t] + GAP),
                                (compare_matrix[k][t - 1] + GAP)
                            );
						}
					}
				}

				#pragma omp critical 
                {
                    scores[scores_count * 3] = compare_matrix[200][200];
					scores[scores_count * 3 + 1] = i;
					scores[scores_count * 3 + 2] = j;
					scores_count += 1;
				}
			}
		}

		#pragma omp barrier

        quick_sort(scores, 0, scores_count * 3 - 3);
		
		for(int i = scores_count * 3 - 3; i > scores_count * 3 - 63; i -=3) {
			max_scores[max_scores_count * 3] = scores[i];
			max_scores[max_scores_count * 3 + 1] = scores[i + 1];
			max_scores[max_scores_count * 3 + 2] = scores[i + 2];
			max_scores_count += 1;
		}
    }

    quick_sort(max_scores, 0, max_scores_count * 3 - 3);

    ofstream result_file(REPORTS);

    result_file << endl << setw(70) << "[SCORE TABLE]" << endl << endl << endl;
    result_file << setw(47) << "No" << setw(10) << "S1" << setw(10) << "S2" << setw(15) << "Score" << endl;
    result_file << setw(86) << "--------------------------------------------" << endl;

    count = 1;
    
    // for(int i = max_scores_count * 3 - 3; i > max_scores_count * 3 - 63; i -=3) {
    //     result_file << setw(47) << count << setw(10) << max_scores[i + 1] << setw(10) << max_scores[i + 2] << setw(17) << max_scores[i] << endl;
    //     result_file << "sequence 1 :" << sequence_arrays[ (int) max_scores[i + 1]] << endl;
    //     result_file << "sequence 2 :" << sequence_arrays[ (int) max_scores[i + 2]] << endl;
    //     result_file << setw(86) << "--------------------------------------------" << endl;
    //     count += 1;
    // }

    result_file.close();

    cout << endl << setw(39) << get_time() << endl;

    return 0;
}

