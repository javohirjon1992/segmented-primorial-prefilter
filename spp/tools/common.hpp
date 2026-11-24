// Small helpers shared by the command-line tools.
#pragma once
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <string>
#include <vector>

namespace tools {
inline double now_ns() {
    return static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::steady_clock::now().time_since_epoch()).count());
}
inline std::string fmt_time(double ns) {
    char b[64];
    if (ns < 1e3) std::snprintf(b, sizeof b, "%.0f ns", ns);
    else if (ns < 1e6) std::snprintf(b, sizeof b, "%.2f us", ns / 1e3);
    else if (ns < 1e9) std::snprintf(b, sizeof b, "%.3f ms", ns / 1e6);
    else std::snprintf(b, sizeof b, "%.3f s", ns / 1e9);
    return b;
}
inline double median(std::vector<double> v) {
    if (v.empty()) return 0;
    std::sort(v.begin(), v.end());
    size_t n = v.size();
    return n % 2 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}
inline std::string shorten(const std::string& s, size_t w = 60) {
    if (s.size() <= w) return s;
    return s.substr(0, 25) + "..." + s.substr(s.size() - 25) + " (" + std::to_string(s.size()) + " xona)";
}
inline std::string group(unsigned long long v) {      // 361328 -> "361 328"
    std::string s = std::to_string(v), out;
    int c = 0;
    for (auto it = s.rbegin(); it != s.rend(); ++it) { if (c && c % 3 == 0) out.insert(out.begin(), ' '); out.insert(out.begin(), *it); ++c; }
    return out;
}
}  // namespace tools
