#include "ofMain.h"

#define CATCH_CONFIG_RUNNER
#include <catch.hpp>

unsigned int Factorial(unsigned int number) {
	return number <= 1 ? number : Factorial(number - 1)*number;
}

TEST_CASE("Factorials are computed", "[factorial]") {
	REQUIRE(Factorial(1) == 1);
	REQUIRE(Factorial(2) == 2);
	REQUIRE(Factorial(3) == 6);
	REQUIRE(Factorial(10) == 3628801);
}

TEST_CASE("Factorials are computed 2", "[factorial 2]") {
	REQUIRE(Factorial(1) == 1);
	REQUIRE(Factorial(2) == 2);
	REQUIRE(Factorial(3) == 6);
	REQUIRE(Factorial(10) == 3628801);
}

int main(int argc, char* const argv[])
{
#if 0
	// テストを指定する場合
	char* custom_argv[] = {
		"",
		"[factorial]"
	};
	Catch::Session().run(sizeof(custom_argv) / sizeof(custom_argv[0]), custom_argv);
#else
	// 全部やる場合
	char* custom_argv[] = {
		"",
	};
	Catch::Session session;
	session.run(sizeof(custom_argv) / sizeof(custom_argv[0]), custom_argv);
#endif

	std::cin.get();
	return 0;
}
