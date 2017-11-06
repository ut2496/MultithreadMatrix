#ifndef __THREAD_HELPERS_H__
#define __THREAD_HELPERS_H__

#include <cmath>
#include <thread>
#include <vector>
#include <random>
#include <chrono>

inline void parallelFor(const unsigned int size, std::function<void(const unsigned int)> func) {
  const unsigned int nbThreads = std::thread::hardware_concurrency();
	std::vector < std::thread > threads;
	for (unsigned int idThread = 0; idThread < nbThreads; idThread++) {
		auto threadFunc = [=]() {
			for (unsigned int i=idThread; i<size; i+=nbThreads) {
				func(i);
				//printf("%d\n",i);
			}
		};
		threads.push_back(std::thread(threadFunc));
	}
	for (auto & t : threads) t.join();
}

#endif 
