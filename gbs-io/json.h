#pragma once
/**
 * @file json.h
 * @brief JSON documents read and written with RapidJSON, without the geometry of fromjson.h and
 * tojson.h (which include this file).
 *
 * Numbers are read exactly: RapidJSON's Writer writes a double with the shortest digits that give
 * it back, and kParseFullPrecisionFlag reads them back bit for bit (the default parser rounds
 * some of them to a neighbour). A syntax error throws JsonError, with where and why. Files are
 * written whole (write_text_atomic).
 */
#include <gbs-io/files.h>

#include <rapidjson/document.h>
#include <rapidjson/error/en.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/writer.h>

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace gbs
{
    /// A JSON text that is not valid: where (offset in bytes from 0, line and column in bytes from
    /// 1) and why (RapidJSON's message). what() gives the source, the place and the reason.
    class JsonError : public std::runtime_error
    {
    public:
        std::size_t offset{}, line{}, column{};
        std::string reason;

        JsonError(std::string_view source, std::size_t offset_, std::size_t line_, std::size_t column_, std::string reason_)
            : std::runtime_error{std::string{source} + ": line " + std::to_string(line_) + ", column " + std::to_string(column_) +
                                 ": " + reason_},
              offset{offset_}, line{line_}, column{column_}, reason{std::move(reason_)}
        {
        }
    };

    /**
     * @brief Parses a JSON text into `document`, numbers in full precision.
     *
     * @param source what the text is, for the message of an error (a file name)
     * @throws JsonError on a syntax error
     */
    inline void parse_json(std::string_view text, rapidjson::Document &document, std::string_view source = "JSON text")
    {
        document.Parse<rapidjson::kParseFullPrecisionFlag>(text.data(), text.size());
        if (!document.HasParseError())
            return;
        const auto offset = std::min(document.GetErrorOffset(), text.size());
        const auto before = text.substr(0, offset);
        const auto line_start = before.rfind('\n');
        throw JsonError{source,
                        offset,
                        static_cast<std::size_t>(std::ranges::count(before, '\n')) + 1,
                        offset - (line_start == std::string_view::npos ? 0 : line_start + 1) + 1,
                        rapidjson::GetParseError_En(document.GetParseError())};
    }

    /**
     * @brief Parses a JSON file into `document`, numbers in full precision.
     *
     * @throws std::runtime_error if the file cannot be read, JsonError on a syntax error
     */
    inline void parse_file(const std::filesystem::path &fname, rapidjson::Document &document)
    {
        const auto text = read_text(fname);
        if (!text)
            throw std::runtime_error("File not found: " + fname.string());
        parse_json(*text, document, fname.string());
    }

    /**
     * @brief Writes a JSON document to a file, compact, whole or not at all (write_text_atomic).
     *
     * @throws std::filesystem::filesystem_error if the file cannot be written
     */
    inline void write_js_doc(const rapidjson::Document &d, const std::filesystem::path &fName)
    {
        rapidjson::StringBuffer buffer;
        rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);
        d.Accept(writer);
        write_text_atomic(fName, {buffer.GetString(), buffer.GetSize()});
    }
}
