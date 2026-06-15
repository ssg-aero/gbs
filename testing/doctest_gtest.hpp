#pragma once
//
// doctest_gtest.hpp — couche de compatibilité GoogleTest -> doctest.
//
// But : permettre de migrer les fichiers de test gbs (qui n'utilisent que
// TEST/TEST_F + un petit jeu d'assertions) vers doctest sans réécrire le corps
// des tests. Il suffit de remplacer `#include <gtest/gtest.h>` par
// `#include <doctest_gtest.hpp>`.
//
// doctest étant header-only, plus de bibliothèque compilée à apparier au
// runtime Debug/Release (fin du couple gtest_dll / gtest_dlld).
//
// Le main n'est PAS défini ici : il est fourni une seule fois par
// testing/doctest_main.cpp, ajouté à chaque exécutable via la lib INTERFACE
// gbs_testing. Cela couvre aussi bien les exés mono-fichier (tests/) que les
// exés mono-binaire multi-fichiers (gbs-occt, tests_3rdlibs).
//
#include <doctest/doctest.h>
#include <memory>

// ---------------------------------------------------------------------------
// Stringification des pointeurs intelligents
// ---------------------------------------------------------------------------
// doctest stringifie les operandes d'une assertion pour le message d'echec.
// Sa detection d'operator<< est SFINAE-friendly, mais l'operator<< que MSVC
// fournit pour shared_ptr/unique_ptr passe la detection puis echoue DUR a
// l'instanciation (`_Out << _Px.get()` : un T* non streamable). On court-circuite
// donc cette voie en stringifiant explicitement par l'adresse pointee, ce qui
// suffit pour comparer des shared_ptr avec ASSERT_EQ/ASSERT_NE.
namespace doctest {
template <typename T>
struct StringMaker<std::shared_ptr<T>> {
    static String convert(const std::shared_ptr<T> &p) {
        return p ? toString(static_cast<const void *>(p.get())) : String("nullptr");
    }
};
template <typename T, typename D>
struct StringMaker<std::unique_ptr<T, D>> {
    static String convert(const std::unique_ptr<T, D> &p) {
        return p ? toString(static_cast<const void *>(p.get())) : String("nullptr");
    }
};
} // namespace doctest

// ---------------------------------------------------------------------------
// Cas de test
// ---------------------------------------------------------------------------
// gtest : TEST(Suite, Name) -> doctest : TEST_CASE("Suite.Name")
#define TEST(suite, name) DOCTEST_TEST_CASE(#suite "." #name)

// Base de fixture facon gtest. SetUp()/TearDown() sont appeles explicitement
// par TEST_F (cf. plus bas), le constructeur/destructeur de la fixture restent
// egalement honores.
namespace testing {
class Test {
public:
    virtual ~Test() = default;
    virtual void SetUp() {}
    virtual void TearDown() {}
};
} // namespace testing

// gtest : TEST_F(Fixture, Name). doctest construit une instance fraiche de la
// fixture par cas ; on encadre le corps par SetUp()/TearDown() pour couvrir les
// fixtures qui s'appuient dessus comme celles a constructeur.
#define TEST_F(fixture, name)                                                   \
    namespace {                                                                 \
    struct GtestF_##fixture##_##name : fixture {                                \
        void gtest_body();                                                      \
    };                                                                          \
    }                                                                           \
    DOCTEST_TEST_CASE(#fixture "." #name) {                                     \
        GtestF_##fixture##_##name _gtest_fixture;                               \
        _gtest_fixture.SetUp();                                                 \
        _gtest_fixture.gtest_body();                                            \
        _gtest_fixture.TearDown();                                              \
    }                                                                           \
    void GtestF_##fixture##_##name::gtest_body()

// ---------------------------------------------------------------------------
// Assertions
// ---------------------------------------------------------------------------
// gtest autorise un message en flux : ASSERT_X(...) << "contexte " << val;
// doctest non. On enveloppe donc chaque assertion pour que le << "..." final
// soit avale par un sink, tout en gardant la macro mono-instruction (idiome
// switch/else de gtest pour neutraliser le "dangling else").
namespace gbs_doctest_compat {
struct MsgSink {
    template <class T> MsgSink &operator<<(const T &) { return *this; }
};

// ASSERT_NEAR : tolerance ABSOLUE. On evite std::abs (ambigu sur les types non
// signes comme uint64_t) en calculant l'ecart sans soustraction signee.
template <class A, class B, class T>
inline bool near(const A &a, const B &b, const T &tol) {
    return (a < b ? b - a : a - b) <= tol;
}
} // namespace gbs_doctest_compat

#define GBS_GTEST_WRAP_(doctest_stmt)                                           \
    switch (0)                                                                  \
    case 0:                                                                     \
    default:                                                                    \
        if (true) {                                                            \
            doctest_stmt;                                                       \
        } else                                                                  \
            ::gbs_doctest_compat::MsgSink {}

// ASSERT_* (fatal) -> REQUIRE* ; EXPECT_* (non fatal) -> CHECK*
// TRUE/FALSE : _UNARY (pas de decomposition) pour accepter les expressions
// booleennes composees (ternaire, &&, ||) que la decomposition doctest refuse.
#define ASSERT_TRUE(cond)  GBS_GTEST_WRAP_(DOCTEST_REQUIRE_UNARY(cond))
#define ASSERT_FALSE(cond) GBS_GTEST_WRAP_(DOCTEST_REQUIRE_UNARY_FALSE(cond))
#define EXPECT_TRUE(cond)  GBS_GTEST_WRAP_(DOCTEST_CHECK_UNARY(cond))
#define EXPECT_FALSE(cond) GBS_GTEST_WRAP_(DOCTEST_CHECK_UNARY_FALSE(cond))

#define ASSERT_EQ(a, b) GBS_GTEST_WRAP_(DOCTEST_REQUIRE((a) == (b)))
#define ASSERT_NE(a, b) GBS_GTEST_WRAP_(DOCTEST_REQUIRE((a) != (b)))
#define ASSERT_LT(a, b) GBS_GTEST_WRAP_(DOCTEST_REQUIRE((a) < (b)))
#define ASSERT_LE(a, b) GBS_GTEST_WRAP_(DOCTEST_REQUIRE((a) <= (b)))
#define ASSERT_GT(a, b) GBS_GTEST_WRAP_(DOCTEST_REQUIRE((a) > (b)))
#define ASSERT_GE(a, b) GBS_GTEST_WRAP_(DOCTEST_REQUIRE((a) >= (b)))

#define EXPECT_EQ(a, b) GBS_GTEST_WRAP_(DOCTEST_CHECK((a) == (b)))
#define EXPECT_NE(a, b) GBS_GTEST_WRAP_(DOCTEST_CHECK((a) != (b)))
#define EXPECT_LT(a, b) GBS_GTEST_WRAP_(DOCTEST_CHECK((a) < (b)))
#define EXPECT_LE(a, b) GBS_GTEST_WRAP_(DOCTEST_CHECK((a) <= (b)))
#define EXPECT_GT(a, b) GBS_GTEST_WRAP_(DOCTEST_CHECK((a) > (b)))
#define EXPECT_GE(a, b) GBS_GTEST_WRAP_(DOCTEST_CHECK((a) >= (b)))

#define ASSERT_NEAR(a, b, tol)                                                  \
    GBS_GTEST_WRAP_(DOCTEST_REQUIRE_UNARY(::gbs_doctest_compat::near((a), (b), (tol))))
#define EXPECT_NEAR(a, b, tol)                                                  \
    GBS_GTEST_WRAP_(DOCTEST_CHECK_UNARY(::gbs_doctest_compat::near((a), (b), (tol))))

// Egalite flottante facon gtest (~4 ULP) -> Approx (epsilon relatif par defaut).
#define ASSERT_DOUBLE_EQ(a, b) GBS_GTEST_WRAP_(DOCTEST_REQUIRE((a) == doctest::Approx(b)))
#define EXPECT_DOUBLE_EQ(a, b) GBS_GTEST_WRAP_(DOCTEST_CHECK((a) == doctest::Approx(b)))
#define ASSERT_FLOAT_EQ(a, b)                                                   \
    GBS_GTEST_WRAP_(DOCTEST_REQUIRE((a) == doctest::Approx(b).epsilon(1e-5)))
#define EXPECT_FLOAT_EQ(a, b)                                                   \
    GBS_GTEST_WRAP_(DOCTEST_CHECK((a) == doctest::Approx(b).epsilon(1e-5)))

#define ASSERT_THROW(stmt, ex) GBS_GTEST_WRAP_(DOCTEST_REQUIRE_THROWS_AS(stmt, ex))
#define EXPECT_THROW(stmt, ex) GBS_GTEST_WRAP_(DOCTEST_CHECK_THROWS_AS(stmt, ex))
#define ASSERT_NO_THROW(stmt)  GBS_GTEST_WRAP_(DOCTEST_REQUIRE_NOTHROW(stmt))
#define EXPECT_NO_THROW(stmt)  GBS_GTEST_WRAP_(DOCTEST_CHECK_NOTHROW(stmt))
#define ASSERT_ANY_THROW(stmt) GBS_GTEST_WRAP_(DOCTEST_REQUIRE_THROWS(stmt))
#define EXPECT_ANY_THROW(stmt) GBS_GTEST_WRAP_(DOCTEST_CHECK_THROWS(stmt))

// ---------------------------------------------------------------------------
// GTEST_SKIP() [<< "message"]
// ---------------------------------------------------------------------------
// doctest 2.4.x n'a pas de skip a chaud au milieu d'un test : on sort du cas
// (il est compte comme reussi, pas comme "skipped"). Le sink avale le << "...".
namespace gbs_doctest_compat {
struct SkipSink {
    template <class T> SkipSink &operator<<(const T &) { return *this; }
};
} // namespace gbs_doctest_compat

#define GTEST_SKIP()                                                            \
    if (gbs_doctest_compat::SkipSink _gtest_skip_sink{}; true)                  \
        return;                                                                 \
    else                                                                        \
        _gtest_skip_sink
