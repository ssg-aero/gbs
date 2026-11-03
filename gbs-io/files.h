#pragma once
/**
 * @file files.h
 * @brief Text files read whole and written whole: no geometry, no third-party library.
 *
 * read_text reads a file in one block (the STEP reader, the JSON files). write_text_atomic
 * writes a file so that a reader never sees it half written: the text goes to a temporary file of
 * a name of its own in the same directory, renamed over the file once complete. utf8_string gives
 * a path for a message on every system.
 */
#include <cerrno>
#include <filesystem>
#include <format>
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

    /// The UTF-8 text of a path, for a message: path::string() converts to the code page of the
    /// system on Windows, and throws for a name it cannot hold.
    inline auto utf8_string(const std::filesystem::path &p) -> std::string
    {
        const auto s = p.u8string();
        return std::string(s.begin(), s.end());
    }

    namespace files_detail
    {
        /// The file a write goes to: `file`, its symbolic links followed (a relative one from its
        /// directory) as an in-place write follows them, at most 40 as the systems do.
        inline auto through_links(std::filesystem::path file) -> std::filesystem::path
        {
            std::error_code e;
            for (int k = 0; k < 40 && std::filesystem::is_symlink(file, e); ++k)
            {
                const auto to = std::filesystem::read_symlink(file, e);
                if (e)
                    break;
                file = file.parent_path() / to; // an absolute target replaces the whole
            }
            return file;
        }

        /// At most n units of a file name, not cut inside a character (a UTF-8 continuation byte,
        /// a UTF-16 low surrogate on Windows)
        inline auto shortened(std::filesystem::path::string_type name, std::size_t n) -> std::filesystem::path::string_type
        {
            if (name.size() <= n)
                return name;
            const auto inside = [](auto c) {
                if constexpr (sizeof(c) == 1)
                    return (static_cast<unsigned char>(c) & 0xC0u) == 0x80u;
                else
                    return c >= 0xDC00 && c <= 0xDFFF;
            };
            while (n > 0 && inside(name[n]))
                --n;
            name.resize(n);
            return name;
        }
    }

    /**
     * @brief Writes a file whole, or leaves it as it was.
     *
     * The text is written to a temporary file next to the file (same directory, so the renaming
     * never copies across file systems), named `.<name>.<16 hex digits>.tmp` from a random draw so
     * that two writers at once never share it (the name cut to 200 units, so that the whole stays
     * within the 255 of a file name), then renamed over the file. A reader opening the file
     * meanwhile gets either the old text or the new one, never a part; of two writers at once, the
     * last renamed wins, whole. On failure the temporary file is removed, the file is untouched,
     * and std::filesystem::filesystem_error is thrown.
     *
     * As an in-place write: a symbolic link is written through (the link stays, its target gets
     * the text), the file keeps its permissions, and a file that cannot be written (read only)
     * is refused. Unlike it, a hard link is broken: the other names of the file keep the old text;
     * the owner and the access control lists are those of the writer.
     *
     * Atomic for the readers, not durable: nothing is flushed to the disk (no fsync), so a power
     * cut just after may lose the new text (not mix it with the old one on journaled file systems).
     *
     * On Windows the renaming (MoveFileExW, MOVEFILE_REPLACE_EXISTING) fails while another process
     * holds the file open without FILE_SHARE_DELETE, as the C and C++ runtimes open files for
     * reading: the call then throws (access denied or sharing violation) and the file keeps its
     * old text. It is not retried: the caller knows whether a reader may hold the file and how
     * long to wait for it.
     */
    inline void write_text_atomic(const std::filesystem::path &file, std::string_view text)
    {
        const auto target = files_detail::through_links(file);
        std::random_device random;
        auto tmp = target.parent_path() / ".";
        tmp += files_detail::shortened(target.filename().native(), 200);
        tmp += std::format(".{:08x}{:08x}.tmp", random(), random());
        // the reason of a failed stream operation: errno, set by the system call under it
        const auto os_error = [] { return std::error_code{errno != 0 ? errno : EIO, std::generic_category()}; };
        const auto fail = [&](const char *what, std::error_code e) {
            std::error_code ignored;
            std::filesystem::remove(tmp, ignored);
            throw std::filesystem::filesystem_error{what, file, tmp, e};
        };
        std::error_code e;
        const auto old = std::filesystem::status(target, e);
        const bool exists = std::filesystem::exists(old);
        // writable, as an in-place write needs it: opened to append, nothing written
        errno = 0;
        if (exists && !std::ofstream{target, std::ios::binary | std::ios::app})
            fail("cannot write the file", os_error());
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
        if (exists)
            std::filesystem::permissions(tmp, old.permissions(), e); // those of the file it replaces
        std::filesystem::rename(tmp, target, e);
        if (e)
            fail("cannot replace the file", e);
    }
}
