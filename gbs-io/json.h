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
#include <cerrno>
#include <cstddef>
#include <filesystem>
#include <format>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>

namespace gbs
{
    /// A JSON text that is not valid: where (offset in bytes from 0, line and column in bytes from
    /// 1) and why (RapidJSON's message, reason()). what() gives the source, the place and the
    /// reason; the reason is the end of it, so that copying a JsonError never throws.
    class JsonError : public std::runtime_error
    {
        std::size_t reason_at_{};

    public:
        std::size_t offset{}, line{}, column{};

        JsonError(std::string_view source, std::size_t offset_, std::size_t line_, std::size_t column_, std::string_view reason_)
            : std::runtime_error{std::format("{}: line {}, column {}: {}", source, line_, column_, reason_)},
              reason_at_{std::string_view{what()}.size() - reason_.size()}, offset{offset_}, line{line_}, column{column_}
        {
        }

        /// why the text is not valid (RapidJSON's message)
        [[nodiscard]] auto reason() const noexcept -> std::string_view { return std::string_view{what()}.substr(reason_at_); }
    };

    /**
     * @brief Parses a JSON text into `document`, numbers in full precision.
     *
     * A NUL byte is an error: RapidJSON ends the text there, so that what follows would not be
     * read. On an error, `document` is null.
     *
     * @param source what the text is, for the message of an error (a file name)
     * @throws JsonError on a syntax error
     */
    inline void parse_json(std::string_view text, rapidjson::Document &document, std::string_view source = "JSON text")
    {
        const auto error_at = [&](std::size_t offset, std::string_view reason) {
            document.SetNull();
            offset = std::min(offset, text.size());
            const auto before = text.substr(0, offset);
            const auto line_start = before.rfind('\n');
            return JsonError{source,
                             offset,
                             static_cast<std::size_t>(std::ranges::count(before, '\n')) + 1,
                             offset - (line_start == std::string_view::npos ? 0 : line_start + 1) + 1,
                             reason};
        };
        document.Parse<rapidjson::kParseFullPrecisionFlag>(text.data(), text.size());
        if (document.HasParseError())
            throw error_at(document.GetErrorOffset(), rapidjson::GetParseError_En(document.GetParseError()));
        if (const auto nul = text.find('\0'); nul != std::string_view::npos)
            throw error_at(nul, "NUL byte in the text.");
    }

    /**
     * @brief Parses a JSON file into `document`, numbers in full precision.
     *
     * @throws std::filesystem::filesystem_error if the file cannot be read (with the reason of the
     * system), JsonError on a syntax error
     */
    inline void parse_file(const std::filesystem::path &fname, rapidjson::Document &document)
    {
        errno = 0;
        const auto text = read_text(fname);
        if (!text)
        {
            const int reason = errno;
            std::error_code e;
            throw std::filesystem::filesystem_error{
                "cannot read the JSON file", fname,
                std::error_code{std::filesystem::is_directory(fname, e) ? EISDIR : reason != 0 ? reason : EIO, std::generic_category()}};
        }
        parse_json(*text, document, utf8_string(fname));
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
