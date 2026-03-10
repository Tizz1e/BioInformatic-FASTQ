#include "fastq-parser.h"

#include <iostream>

bool operator==(const Record& lhs, const Record& rhs) {
    return std::tie(lhs.label, lhs.raw, lhs.description, lhs.quality) == std::tie(rhs.label, rhs.raw, rhs.description, rhs.quality);
}

void Print(const Record& rec) {
    std::cout << rec.label << '\n';
    std::cout << rec.raw << '\n';
    std::cout << rec.description << '\n';
    std::cout << rec.quality << '\n';
}

FastqParser::FastqParser(const std::filesystem::path& in) : in_(in) {
    GetNext();
}

bool FastqParser::HasNextRecord() const {
    return record_.has_value();
}

Record FastqParser::NextRecord() {
    auto rec = std::move(*record_);
    GetNext();
    return rec;
}

void FastqParser::GetNext() {
    if (!in_.is_open()) {
        return;
    }
    record_.reset();
    
    Record rec;
    std::string line;

    getline(in_, line);
    rec.label = std::move(line);

    getline(in_, line);
    rec.raw = std::move(line);

    getline(in_, line);
    rec.description = std::move(line);

    getline(in_, line);
    rec.quality = std::move(line);

    if (in_.eof()) {
        in_.close();
        return;
    }
    
    record_ = rec;
}
