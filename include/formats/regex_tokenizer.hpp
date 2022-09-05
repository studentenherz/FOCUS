#if !defined(FOCUS_INCLUDE_FORMATS_REGEX_TOKENIZER_HPP)
#define FOCUS_INCLUDE_FORMATS_REGEX_TOKENIZER_HPP

#include <regex>
#include <string>

/**
 * Class to tokenize a stream based on a regex pattern
 */
template <typename stream_type>
class Tokenizer{
	std::regex e;
	std::string line;
	std::smatch m;
public:
	/**
	 * @param rexp Regular expresion
	 */
	Tokenizer(const char rexp[]): e(rexp) {}

	/**
	 * Get next token from stream
	 * @param is stream to tokenize
	 * @param[out] token string to put the token in
	 * @return true if new token, false otherwise
	 */
	bool next(stream_type& is, std::string& token){
		if (!std::regex_search (line,m,e)){
			if (!std::getline(is, line))
				return false;
			if (!std::regex_search (line,m,e))
				return false;
		}
    token = m[0];
		line = m.suffix().str();
		return true;
	}
};

/**
 * Check if string matches a regular expresion
 * @param s string
 * @param rexpr regular expresion
 * @return true if s matches rexpr
 */
inline bool regex_match(std::string s, std::string rexpr){
	std::regex e(rexpr);
	return std::regex_match(s, e);
}

#endif // FOCUS_INCLUDE_FORMATS_REGEX_TOKENIZER_HPP
