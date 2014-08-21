/*
 * Perf.cpp
 *
 *  Created on: 13 févr. 2014
 *      Author: Adrien CASSAGNE
 */

#include <sys/time.h>
#include <cassert>
#include <iostream>

#include "Perf.h"

Perf::Perf()
	: tStart(0), tStop(0) {
}

Perf::Perf(const Perf &p)
	: tStart(p.tStart), tStop(p.tStop) {
}

Perf::Perf(double ms)
	: tStart(0), tStop(ms * 1000) {
}

Perf::~Perf() {
}

void Perf::start() {
    this->tStart = Perf::getTime();
}

void Perf::stop() {
	this->tStop = Perf::getTime();
}

void Perf::reset() {
	this->tStart = 0;
	this->tStop  = 0;
}

double Perf::getElapsedTime() {
	return (this->tStop - this->tStart) / 1000.f;
}

double Perf::getGflops(unsigned long flops) {
	return (flops * (1000 / this->getElapsedTime())) / 1024.0 / 1024.0 / 1024.0;
}

double Perf::getMemoryBandwidth(unsigned long memops, unsigned short nBytes) {
	return (memops * nBytes * (1000 / this->getElapsedTime())) / 1024.0 / 1024.0 / 1024.0;
}

Perf Perf::operator+(const Perf& p) {
	Perf pAdd;
	pAdd.tStop = (p.tStop - p.tStart) + (this->tStop - this->tStart);
	return pAdd;
}

Perf Perf::operator+=(const Perf& p) {
	this->tStop += p.tStop - p.tStart;
	return (*this);
}

unsigned long Perf::getTime() {
    struct timeval t;

    int ret = gettimeofday(&t, NULL);
    assert(ret == 0);

    return t.tv_sec * 1000000 + t.tv_usec;
}
