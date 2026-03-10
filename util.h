#pragma once

#include "fastq-parser.h"
#include <limits>

struct Stat {
    size_t total_reads = 0;
    size_t min_read = std::numeric_limits<size_t>::max();
    size_t max_read = 0;
    size_t avg_read = 0;
};

Stat GetStatistics(const std::filesystem::path& path);
double GCContent(const std::filesystem::path& path);
void Trim(const std::filesystem::path& in, const std::filesystem::path& out, std::size_t width, double quality);
void Print(const Stat& stat);