/*
 * Perf.h
 *
 *  Created on: 13 févr. 2014
 *      Author: Adrien CASSAGNE
 */

#ifndef PERF_H_
#define PERF_H_

class Perf {
private:
	unsigned long tStart;
	unsigned long tStop;

public:
	Perf();
	Perf(const Perf &p);
	Perf(float ms);
	virtual ~Perf();

	void start();
	void stop();
	void reset();

	float getElapsedTime();                                                // ms
	float getGflops(float flops);                                          // Gflops/s
	float getMemoryBandwidth(unsigned long memops, unsigned short nBytes); // Go/s

	Perf operator+(const Perf& p);
	Perf operator+=(const Perf& p);

protected:
	static unsigned long getTime();
};

#endif /* PERF_H_ */
