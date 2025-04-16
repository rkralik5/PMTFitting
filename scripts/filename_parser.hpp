#pragma once

#include <string>
#include <vector>
#include <optional>
#include <stdexcept>  // for std::invalid_argument

///
/// Represent PMT information from the measurement filename.
/// If multiple channels are used, the `channel` will be set (e.g. "Ch1").
///
struct FileInfo {
  std::optional<float> channel;
  std::string manufacturer;
  std::string model;
  std::string vstring;
  float voltage;
};

/// Holds the result of parsing a PMT measurement filename.
struct ParsedFile {
  std::string path;
  std::string filename;
  std::vector<FileInfo> devices;
  std::string description;
};

///
/// Splits a string into tokens based on a specified delimiter.
///
/// @param input The input string to be split
/// @param delimiter The character used to split the string (default is '_')
/// @return A vector of strings containing the tokens
///
std::vector<std::string> splitString(const std::string& input,
                                     char delimiter = '_') {
  std::vector<std::string> tokens;
  std::string current;
  for (char c : input) {
    if (c == delimiter) {
      if (!current.empty()) tokens.push_back(current);
      current.clear();
    } else {
      current += c;
    }
  }
  if (!current.empty()) tokens.push_back(current);
  return tokens;
}

///
/// Parses a PMT measurement filename and extracts structured information.
/// 
/// @param filename The full filename, e.g. "Ch1_Hamamatsu_R12860_1300V.root"
/// @return ParsedFile object containing extracted metadata
///
ParsedFile parseFilename(const std::string& fullFilename) {
  ParsedFile result;

  // Extract the filename without the path
  result.path = fullFilename.substr(0, fullFilename.find_last_of('/'));
  result.filename = fullFilename.substr(fullFilename.find_last_of("/")+1);
  std::string baseName = result.filename.substr(0, result.filename.find_last_of('.'));

  // Split the base name into tokens using underscore as a delimiter
  std::vector<std::string> tokens = splitString(baseName);

  // Loop over channels and extract information
  // The regex checks if the first token is a channel (e.g. "Ch1")
  // If it is, we expect the next three tokens to be manufacturer, model, and
  // voltage, otherwise they are the first three tokens
  size_t i = 0;
  while (i + 3 <= tokens.size()) {
    FileInfo info;

    if (tokens[i].substr(0, 2) == "Ch") {
      info.channel = std::stof(tokens[i].substr(2));
      info.manufacturer = tokens[i+1];
      info.model = tokens[i+2];
      info.vstring = tokens[i+3];
      if (std::tolower(tokens[i+3].back()) != 'v') {
        throw std::invalid_argument("Invalid voltage format: "+tokens[i+3]);
      }
      info.voltage = std::stof(tokens[i+3].substr(0, tokens[i+3].size()-1));
      i += 4;
    } else {
      info.manufacturer = tokens[i];
      info.model = tokens[i+1];
      info.vstring = tokens[i+2];
      if (std::tolower(tokens[i+2].back()) != 'v') {
        throw std::invalid_argument("Invalid voltage format: "+tokens[i+2]);
      }
      info.voltage = std::stof(tokens[i+2].substr(0, tokens[i+2].size()-1));
      i += 3;
    }

    result.devices.push_back(info);
  }

  // If there are any remaining tokens, they are part of the description
  if (i < tokens.size()) {
    std::string desc = tokens[i];
    for (size_t j = i + 1; j < tokens.size(); ++j) {
      desc += "_" + tokens[j];
    }
    result.description = desc;
  } else {
    result.description = "";
  }

  return result;
}