#include <iostream>
#include <algorithm>
#include <random>
#include <chrono>
#include <iomanip>
#include <thread>
#include <mutex>
using namespace std;

std::mutex consoleMutex;

int barrierSearch(int* arr, int size, int key) {
    int last = arr[size - 1];
    arr[size - 1] = key;

    int i = 0;
    while (arr[i] != key) {
        i++;
    }

    arr[size - 1] = last;

    if (i < size - 1 || key == last) {
        return i;
    } else {
        return -1;
    }
}

int interpolationSearch(int* arr, int size, int key) {
    int low = 0;
    int high = size - 1;

    while (low <= high && key >= arr[low] && key <= arr[high]) {
        int pos = low + ((static_cast<double>(key - arr[low]) * (high - low)) /
                         (arr[high] - arr[low]));

        if (pos < 0 || pos >= size)
            return -1;

        if (arr[pos] == key) {
            return pos;
        }

        if (arr[pos] < key) {
            low = pos + 1;
        } else {
            high = pos - 1;
        }
    }

    return -1;
}

template <typename Func>
double measureTime(Func func, int iterations) {
    auto start = chrono::high_resolution_clock::now();

    for (int i = 0; i < iterations; i++) {
        func();
    }

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double, milli> duration = end - start;

    return duration.count() / iterations;
}

void parallelSearch(int* arr, int size, int* keys, int keyCount, int* results, int threadId, int numThreads, bool isBarrier) {
    int keysPerThread = keyCount / numThreads;
    int startIdx = threadId * keysPerThread;
    int endIdx = (threadId == numThreads - 1) ? keyCount : (threadId + 1) * keysPerThread;

    for (int i = startIdx; i < endIdx; i++) {
        if (isBarrier) {
            results[i] = barrierSearch(arr, size, keys[i]);
        } else {
            results[i] = interpolationSearch(arr, size, keys[i]);
        }
    }
}

void runExperiment(int size) {
    random_device rd;
    mt19937 gen(rd());

    int* data = new int[size];
    uniform_int_distribution<> dis(1, size * 10);
    for (int i = 0; i < size; i++) {
        data[i] = dis(gen);
    }
    sort(data, data + size);

    const int SEARCH_COUNT = 100;

    int* searchKeys = new int[SEARCH_COUNT];
    for (int i = 0; i < SEARCH_COUNT; i++) {
        int idx = uniform_int_distribution<>(0, size - 1)(gen);
        searchKeys[i] = data[idx];
    }

    int* notFoundKeys = new int[SEARCH_COUNT];
    for (int i = 0; i < SEARCH_COUNT; i++) {
        notFoundKeys[i] = size * 10 + i;
    }

    int* barrierData = new int[size];
    copy(data, data + size, barrierData);

    int* foundResults = new int[SEARCH_COUNT];
    int* notFoundResults = new int[SEARCH_COUNT];

    unsigned int numThreads = min(4u, thread::hardware_concurrency());
    if (numThreads == 0) numThreads = 2;

    thread* barrierFoundThreads = new thread[numThreads];
    thread* barrierNotFoundThreads = new thread[numThreads];
    thread* interpolationFoundThreads = new thread[numThreads];
    thread* interpolationNotFoundThreads = new thread[numThreads];

    double barrierTimeFound = measureTime([&]() {
        for (unsigned int t = 0; t < numThreads; t++) {
            barrierFoundThreads[t] = thread(parallelSearch, barrierData, size,
                                            searchKeys, SEARCH_COUNT, foundResults,
                                            t, numThreads, true);
        }

        for (unsigned int t = 0; t < numThreads; t++) {
            barrierFoundThreads[t].join();
        }
    }, 10);

    double barrierTimeNotFound = measureTime([&]() {
        for (unsigned int t = 0; t < numThreads; t++) {
            barrierNotFoundThreads[t] = thread(parallelSearch, barrierData, size,
                                              notFoundKeys, SEARCH_COUNT, notFoundResults,
                                              t, numThreads, true);
        }

        for (unsigned int t = 0; t < numThreads; t++) {
            barrierNotFoundThreads[t].join();
        }
    }, 10);

    double interpolationTimeFound = measureTime([&]() {
        for (unsigned int t = 0; t < numThreads; t++) {
            interpolationFoundThreads[t] = thread(parallelSearch, data, size,
                                                 searchKeys, SEARCH_COUNT, foundResults,
                                                 t, numThreads, false);
        }

        for (unsigned int t = 0; t < numThreads; t++) {
            interpolationFoundThreads[t].join();
        }
    }, 10);

    double interpolationTimeNotFound = measureTime([&]() {
        for (unsigned int t = 0; t < numThreads; t++) {
            interpolationNotFoundThreads[t] = thread(parallelSearch, data, size,
                                                    notFoundKeys, SEARCH_COUNT, notFoundResults,
                                                    t, numThreads, false);
        }

        for (unsigned int t = 0; t < numThreads; t++) {
            interpolationNotFoundThreads[t].join();
        }
    }, 10);

    {
        lock_guard<mutex> lock(consoleMutex);
        cout << "Array size: " << size << " (threads used: " << numThreads << ")" << endl;
        cout << "----------------------------------------\n";
        cout << "| Search Method       | Found (ms) | Not Found (ms) |\n";
        cout << "----------------------------------------\n";
        cout << "| Barrier Search      | " << setw(11) << fixed << setprecision(4) << barrierTimeFound
              << " | " << setw(14) << fixed << setprecision(4) << barrierTimeNotFound << " |\n";
        cout << "| Interpolation       | " << setw(11) << fixed << setprecision(4) << interpolationTimeFound
              << " | " << setw(14) << fixed << setprecision(4) << interpolationTimeNotFound << " |\n";
        cout << "----------------------------------------\n\n";
    }

    delete[] data;
    delete[] barrierData;
    delete[] searchKeys;
    delete[] notFoundKeys;
    delete[] foundResults;
    delete[] notFoundResults;
    delete[] barrierFoundThreads;
    delete[] barrierNotFoundThreads;
    delete[] interpolationFoundThreads;
    delete[] interpolationNotFoundThreads;
}

int main() {
    cout << "Comparing search methods: Barrier Search vs Interpolation Search\n";
    cout << "Using multiple threads for faster search\n\n";

    thread t1(runExperiment, 1000);
    thread t2(runExperiment, 10000);
    thread t3(runExperiment, 100000);

    t1.join();
    t2.join();
    t3.join();

    return 0;
}

} // end namespace std