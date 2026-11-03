#pragma once
/**
 * @file files.h
 * @brief Text files read whole and written whole: no geometry, no third-party library.
 *
 * read_text reads a file in one block (the STEP reader, the JSON files). write_text_atomic
 * writes a file so that a reader never sees it half written: the text goes to a temporary file of
 * a name of its own in the same directory, renamed over the file once complete.
 */
#include <cerrno>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <optional>
#include <random>
#include <string>
#include <string_view>
#include <system_error>

namespace gbs
{
    /// The whole text of a file, read in one block (a pipe or a special file read as a stream);
    /// none if it cannot be opened or read, or is a directory.
    inline auto read_text(const std::filesystem::path &file) -> std::optional<std::string>
    {
        if (std::error_code ec; std::filesystem::is_directory(file, ec))
            return std::nullopt;
        std::ifstream in(file, std::ios::binary);
        if (!in)
            return std::nullopt;
        std::string text;
        in.seekg(0, std::ios::end);
        if (const auto end = in.tellg(); end != std::streampos(-1) && in.seekg(0, std::ios::beg))
        {
            text.resize(static_cast<std::size_t>(end));
            in.read(text.data(), static_cast<std::streamsize>(text.size()));
            text.resize(static_cast<std::size_t>(in.gcount())); // in case the file shrank meanwhile
        }
        else // pipe or special file
        {
            in.clear();
            text.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
        }
        if (in.bad())
            return std::nullopt;
        return text;
    }

    /**
     * @brief Writes a file whole, or leaves it as it was.
     *
     * The text is written to a temporary file next to `file` (same directory, so the renaming
     * never copies across file systems), named `.<name>.<16 hex digits>.tmp` from a random draw so
     * that two writers at once never share it, then renamed over `file`. A reader opening `file`
     * meanwhile gets either the old text or the new one, never a part; of two writers at once, the
     * last renamed wins, whole. On failure the temporary file is removed, `file` is untouched, and
     * std::filesystem::filesystem_error is thrown.
     *
     * Atomic for the readers, not durable: nothing is flushed to the disk (no fsync), so a power
     * cut just after may lose the new text (not mix it with the old one on journaled file systems).
     *
     * On Windows the renaming (MoveFileExW, MOVEFILE_REPLACE_EXISTING) fails while another process
     * holds `file` open without FILE_SHARE_DELETE, as the C and C++ runtimes open files for
     * reading: the call then throws (access denied or sharing violation) and `file` keeps its old
     * text. It is not retried: the caller knows whether a reader may hold the file and how long to
     * wait for it.
     */
    inline void write_text_atomic(const std::filesystem::path &file, std::string_view text)
    {
        std::random_device random;
        auto tmp = file.parent_path() / ("." + file.filename().string() + ".");
        constexpr char digits[] = "0123456789abcdef";
        for (int k = 0; k < 2; ++k)
            for (auto r = random(), i = 0u; i < 8; ++i, r >>= 4)
                tmp += digits[r & 0xf];
        tmp += ".tmp";
        // the reason of a failed stream operation: errno, set by the system call under it
        const auto os_error = [] { return std::error_code{errno != 0 ? errno : EIO, std::generic_category()}; };
        const auto fail = [&](const char *what, std::error_code e) {
            std::error_code ignored;
            std::filesystem::remove(tmp, ignored);
            throw std::filesystem::filesystem_error{what, file, tmp, e};
        };
        {
            errno = 0;
            std::ofstream out(tmp, std::ios::binary | std::ios::trunc);
            if (!out)
                fail("cannot create the temporary file", os_error());
            out.write(text.data(), static_cast<std::streamsize>(text.size()));
            out.close();
            if (!out)
                fail("cannot write the temporary file", os_error());
        }
        std::error_code e;
        std::filesystem::rename(tmp, file, e);
        if (e)
            fail("cannot replace the file", e);
    }
}
