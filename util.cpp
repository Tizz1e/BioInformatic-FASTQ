#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <optional>

#include "fastq-parser.h"
#include "util.h"

#include <cassert>

namespace {
std::optional<Record>
TrimRecord(Record rec, std::size_t width, double quality) {
  if (rec.quality.size() < width)
    return std::nullopt;

  int total = 0;
  for (int i = 0; i < width; i++) {
    assert(rec.quality[i] - 33 >= 0);
    total += rec.quality[i] - 33;
  }

  if (static_cast<double>(total) < quality * width)
    return std::nullopt;

  int lengthToKeep = rec.quality.size();

  for (int i = 0; i < rec.quality.size() - width; i++) {
    total = total - (rec.quality[i] - 33) + (rec.quality[i + width] - 33);
    if (static_cast<double>(total) < quality * width) {
      lengthToKeep = i + width;
      break;
    }
  }

  int i = lengthToKeep;

  int lastBaseQuality = rec.quality[i - 1] - 33;
  while (lastBaseQuality < quality && i > 1) {
    i--;
    lastBaseQuality = rec.quality[i - 1] - 33;
  }

  if (i < 1)
    return std::nullopt;

  if (i < rec.quality.size()) {
    rec.raw = rec.raw.substr(0, i);
    rec.quality = rec.quality.substr(0, i);
  }

  return rec;
}
} // namespace

Stat GetStatistics(const std::filesystem::path &path) {
  FastqParser parser(path);

  Stat statistics;

  std::size_t sum_length = 0;
  while (parser.HasNextRecord()) {
    auto rec = parser.NextRecord();
    statistics.total_reads++;
    statistics.min_read = std::min(statistics.min_read, rec.raw.size());
    statistics.max_read = std::max(statistics.max_read, rec.raw.size());
    sum_length += rec.raw.size();
  }

  if (statistics.total_reads == 0) {
    return statistics;
  }

  statistics.avg_read = sum_length / statistics.total_reads;

  return statistics;
}

double GCContent(const std::filesystem::path &path) {
  FastqParser parser(path);

  std::size_t gc_count = 0;
  std::size_t total_count = 0;

  while (parser.HasNextRecord()) {
    auto rec = parser.NextRecord();
    gc_count += std::count_if(rec.raw.begin(), rec.raw.end(), [](char c) {
      c = std::toupper(c);
      return c == 'G' || c == 'C';
    });
    total_count += rec.raw.size();
  }

  return 100 * static_cast<double>(gc_count) / total_count;
}

void Trim(const std::filesystem::path &in,
                 const std::filesystem::path &out, std::size_t width,
                 double quality) {
  FastqParser parser(in);
  std::ofstream out_file(out);

  // std::cout << "IN\n";
  while (parser.HasNextRecord()) {
    // std::cout << "PARSING\n";
    auto rec = parser.NextRecord();
    // Print(rec);
    // if (find(rec.raw.begin(), rec.raw.end(), 'N') != rec.raw.end()) {
    //   // std::cout << "FOUND N IN GEN\n";
    //   continue;
    // }


    auto opt_trimmed = TrimRecord(rec, width, quality);
    if (opt_trimmed) {
      auto& trimmed = *opt_trimmed;
      out_file << trimmed.label << '\n';
      out_file << trimmed.raw << '\n';
      out_file << trimmed.description << '\n';
      out_file << trimmed.quality << '\n';
      // std::cout << "GSDSDG" << std::endl;
    } else {
      // std::cout << "NULL\n" << std::endl;
    }
  }
}

void Print(const Stat &stat) {
  std::cout << "Total records: " << stat.total_reads << '\n';
  std::cout << "Min length: " << stat.min_read << '\n';
  std::cout << "Max length: " << stat.max_read << '\n';
  std::cout << "Avg length: " << stat.avg_read << '\n';
}