#pragma once

#include <filesystem>
#include <fstream>
#include <optional>

struct Record {
    std::string label;
    std::string raw;
    std::string description;
    std::string quality;
};

bool operator==(const Record& lhs, const Record& rhs);

void Print(const Record& rec);

class FastqParser {
public:
    FastqParser(const std::filesystem::path& in);
    
    bool HasNextRecord() const;
    Record NextRecord();
private:
    void GetNext();

    std::ifstream in_;
    std::optional<Record> record_;
};