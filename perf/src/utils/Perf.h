/*
 * Do not remove.
 * Optimization training courses 2014 (CINES)
 * Adrien Cassagne, adrien.cassagne@cines.fr
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
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
