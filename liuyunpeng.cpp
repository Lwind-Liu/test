#include <bits/stdc++.h>

#define int long long

using namespace std;

using i64 = long long;

using _T = double;

// constexpr _T eps = 0;
constexpr _T eps = 1e-8;
constexpr long double PI = 3.1415926535897932384l;

template<typename T>
struct point {
	T x, y;

	point(T _x = 0, T _y = 0) : x(_x), y(_y) {}

	bool operator==(const point &p) const { return abs(x - p.x) <= eps && abs(y - p.y) <= eps; }
	bool operator!=(const point &p) const { return !((*this) == p); }
	bool operator<(const point &p) const { return abs(x - p.x) <= eps ? y < p.y - eps : x < p.x - eps; }
	bool operator>(const point &p) const { return !(*this < p || *this == p); }
	point operator+(const point &p) const { return point(x + p.x, y + p.y); }
	point operator-(const point &p) const { return point(x - p.x, y - p.y); }
	point operator-() const { return point(-x, -y); }
	point operator*(const T k) const { return point(x * k, y * k); }
	point operator/(const T k) const { return point(x / k, y / k); }
	T operator*(const point &p) const { return x * p.x + y * p.y; } // Dot
	T operator^(const point &p) const { return x * p.y - y * p.x; } // Cross
	// 1: left  -1: right  0: in
	int toleft(const point &p) const { const auto t = (*this) ^ p; return (t > eps) - (t < -eps); }
	T len2() const { return (*this) * (*this); }
	T dis2(const point &p) const { return (p - (*this)).len2(); }

    // double 相关
	double len() const { return sqrt(len2()); }
	double dis(const point &p) const { return sqrt(dis2(p)); }
    // 夹角
	double ang(const point &p) const { return acos(max(-1.0, min(1.0, ((*this) * p) / len() * p.len()))); }
	// 逆时针旋转rad
    point rot(const double rad) const { return point(x * cos(rad) - y * sin(rad), x * sin(rad) + y * cos(rad)); }
	point rot(const double cosr, const double sinr) const { return point(x * cosr - y * sinr, x * sinr + y * cosr); }
    friend istream& operator>> (istream& is, point& p) { cin >> p.x >> p.y; return is; }
    friend ostream& operator<< (ostream& out, const point& p) { out << p.x << " " << p.y; return out; }
};

using Point = point<_T>;
using Vector = point<_T>;

// 极角序
struct argcmp {
    bool operator()(const Point &a,const Point &b) const {
        const auto quad=[](const Point &a) {
            if (a.y < -eps) return 1;
            if (a.y > eps) return 4;
            if (a.x < -eps) return 5;
            if (a.x > eps) return 3;
            return 2;
        };
        const int qa = quad(a), qb = quad(b);
        if (qa != qb) return qa < qb;
        const auto t = a ^ b;
        // if (abs(t) <= eps) return a * a < b * b - eps;  // 不同长度的向量需要分开
        return t > eps;
    }
};

template<typename T>
struct line {
	point<T> p, v; // p 为一点，v 为方向向量

	line(point<T> _p = point<T>(), point<T> _v = point<T>()) : p(_p), v(_v) {}

	bool operator==(const line &a) const { return v.toleft(a.v) == 0 && v.toleft(p - a.p) == 0; }
	int toleft(const point<T> &a) const { return v.toleft(a - p); }
    bool operator<(const line &a) const { // 半平面交算法定义的排序
        if (abs(v ^ a.v) <= eps && v * a.v >= -eps) return toleft(a.p) == -1;
        return argcmp()(v, a.v);
    }

	// 直线交点
	point<T> inter(const line &a) const { return p + v * ((a.v ^ (p - a.p)) / (v ^ a.v)); }
	// 点到直线距离
    double dis(const point<T> &a) const { return abs(v ^ (a - p)) / v.len(); }
	// a在直线上投影
    point<T> proj(const point<T> &a) const { return p + v * ((v * (a - p)) / (v * v)); }
};

using Line = line<_T>;

// 线段
template<typename T> 
struct segment {
    point<T> a, b;

    segment(point<T> _a = point<T>(), point<T> _b = point<T>()) : a(_a), b(_b) {}
    
    bool operator<(const segment &s) const { return make_pair(a, b) < make_pair(s.a, s.b); }

    // 判定性函数建议在整数域使用
    // 判断点是否在线段上
    // -1 点在线段端点 | 0 点不在线段上 | 1 点严格在线段上
    int is_on(const point<T> &p) const   {
        if (p == a || p == b) return -1;
        return (p - a).toleft(p - b) == 0 && (p - a) * (p - b) < -eps;
    }

    // 判断线段直线是否相交
    // -1 直线经过线段端点 | 0 线段和直线不相交 | 1 线段和直线严格相交
    int is_inter(const line<T> &l) const {
        if (l.toleft(a) == 0 || l.toleft(b) == 0) return -1;
        return l.toleft(a) != l.toleft(b);
    }
    
    // 判断两线段是否相交
    // -1 在某一线段端点处相交 | 0 两线段不相交 | 1 两线段严格相交
    int is_inter(const segment<T> &s) const {
        if (is_on(s.a) || is_on(s.b) || s.is_on(a) || s.is_on(b)) return -1;
        const line<T> l(a, b - a), ls(s.a, s.b - s.a);
        return l.toleft(s.a) * l.toleft(s.b) == -1 && ls.toleft(a) * ls.toleft(b) == -1;
    }

    // 点到线段距离
    long double dis(const point<T> &p) const {
        if ((p - a) * (b - a) < -eps || (p - b) * (a - b) < -eps) return min(p.dis(a), p.dis(b));
        const line<T> l(a, b - a);
        return l.dis(p);
    }

    // 两线段间距离
    long double dis(const segment<T> &s) const {
        if (is_inter(s)) return 0;
        return min({dis(s.a), dis(s.b), s.dis(a), s.dis(b)});
    }

    line<T> to_line() { return line<T>(a, b - a); }
};

using Segment = segment<_T>;

// 多边形
template <typename T>
struct polygon {
    vector<point<T>> p; // 以逆时针顺序存储

    polygon(int sz = 0) : p(sz) {};
    polygon(vector<point<T>> &_p) : p(_p) {}

    size_t nxt(const size_t i) const { return i == p.size() - 1 ? 0 : i + 1; }
    size_t pre(const size_t i) const { return i == 0 ? p.size() - 1 : i - 1; }

    // 回转数
    // 返回值第一项表示点是否在多边形边上
    // 对于狭义多边形，回转数为 0 表示点在多边形外，否则点在多边形内
    pair<bool, int> winding(const point<T> &a) const {
        int cnt = 0;
        for (size_t i = 0; i < p.size(); i ++ ) {
            const point<T> u = p[i], v = p[nxt(i)];
            if (abs((a - u) ^ (a - v)) <= eps && (a - u) * (a - v) <= eps) { return {true, 0}; }
            if (abs(u.y - v.y) <= eps) { continue; }
            const Line uv = {u, v - u};
            if (u.y < v.y - eps && uv.toleft(a) <= 0) { continue; }
            if (u.y > v.y + eps && uv.toleft(a) >= 0) { continue; } 
            if (u.y < a.y - eps && v.y >= a.y - eps) { cnt ++ ; }
            if (u.y >= a.y - eps && v.y < a.y - eps) { cnt -- ; }
        }
        return {false, cnt};
    }

    // 多边形面积的两倍
    // 可用于判断点的存储顺序是顺时针或逆时针
    T area() const {
        T sum = 0;
        for (size_t i = 0; i < p.size(); i ++ ) {
            sum += p[i] ^ p[nxt(i)];
        }
        return sum;
    }

    // 多边形的周长
    long double circ() const {
        long double sum = 0;
        for (size_t i = 0; i < p.size(); i ++ ) {
            sum += p[i].dis(p[nxt(i)]);
        }
        return sum;
    }

    point<T>& operator[](int x) { return p[x]; }
};

using Polygon = polygon<_T>;

// 凸多边形
template <typename T>
struct convex : polygon<T> {
    convex(int sz = 0) : polygon<T>(sz) {}
    convex(vector<point<T>> &_p) : polygon<T>(_p) {}

    // 闵可夫斯基和
    convex operator+(const convex &c) const {
        const auto &p = this->p;
        vector<Segment> e1(p.size()), e2(c.p.size()), edge(p.size() + c.p.size());
        vector<point<T>> res;
        res.reserve(p.size() + c.p.size());
        const auto cmp = [](const Segment &u, const Segment &v) { return argcmp()(u.b - u.a, v.b - v.a); };
        for (size_t i = 0; i < p.size(); i ++ ) { e1[i] = {p[i], p[this->nxt(i)]}; }
        for (size_t i = 0; i < c.p.size(); i++) { e2[i] = {c.p[i], c.p[c.nxt(i)]}; }
        rotate(e1.begin(), min_element(e1.begin(), e1.end(), cmp), e1.end());
        rotate(e2.begin(), min_element(e2.begin(), e2.end(), cmp), e2.end());
        merge(e1.begin(), e1.end(), e2.begin(), e2.end(), edge.begin(), cmp);
        const auto check = [](const vector<point<T>> &res, const point<T> &u) {
            const auto back1 = res.back(), back2 = *prev(res.end(), 2);
            return (back1 - back2).toleft(u - back1) == 0 && (back1 - back2) * (u - back1) >= -eps;
        };
        auto u = e1[0].a + e2[0].a;
        for (const auto &v : edge) {
            while (res.size() > 1 && check(res, u)) {
                res.pop_back();
            }
            res.push_back(u);
            u = u + v.b - v.a;
        }
        if (res.size() > 1 && check(res, res[0])) {
            res.pop_back();
        }
        return { res };
    }

    // 旋转卡壳
    // func 为更新答案的函数，可以根据题目调整位置
    template <typename F>
    void rotcaliper(const F &func) const {
        const auto &p = this->p;
        const auto area = [](const point<T> &u, const point<T> &v, const point<T> &w) { 
            return (w - u) ^ (w - v); 
        };
        for (size_t i = 0, j = 1; i < p.size(); i ++ ) {
            const auto nxti = this->nxt(i);
            func(p[i], p[nxti], p[j]);
            while (area(p[this->nxt(j)], p[i], p[nxti]) >= area(p[j], p[i], p[nxti])) {
                j = this->nxt(j);
                func(p[i], p[nxti], p[j]);
            }
        }
    }

    // 凸多边形的直径的平方
    T diameter2() const {
        const auto &p = this->p;
        if (p.size() == 1) { return 0; }
        if (p.size() == 2) { return p[0].dis2(p[1]); }
        T ans = 0;
        auto func = [&](const point<T> &u, const point<T> &v, const point<T> &w) { 
            ans = max({ans, w.dis2(u), w.dis2(v)}); 
        };
        rotcaliper(func);
        return ans;
    }

    // 判断点是否在凸多边形内
    // 复杂度 O(logn)
    // -1 点在多边形边上 | 0 点在多边形外 | 1 点在多边形内
    int is_in(const point<T> &a) const {
        const auto &p = this->p;
        if (p.size() == 1) { return a == p[0] ? -1 : 0; }
        if (p.size() == 2) { return segment<T>{p[0], p[1]}.is_on(a) ? -1 : 0; }
        if (a == p[0]) { return -1; }
        if ((p[1] - p[0]).toleft(a - p[0]) == -1 || (p.back() - p[0]).toleft(a - p[0]) == 1) { return 0; }
        const auto cmp = [&](const Point &u, const Point &v) { 
            return (u - p[0]).toleft(v - p[0]) == 1; 
        };
        const size_t i = lower_bound(p.begin() + 1, p.end(), a, cmp) - p.begin();
        if (i == 1) { return segment<T>{p[0], p[i]}.is_on(a) ? -1 : 0; }
        if (i == p.size() - 1 && segment<T>{p[0], p[i]}.is_on(a)) { return -1; }
        if (segment<T>{p[i - 1], p[i]}.is_on(a)) { return -1; }
        return (p[i] - p[i - 1]).toleft(a - p[i - 1]) > 0;
    }

    // 凸多边形关于某一方向的极点
    // 复杂度 O(logn)
    // 参考资料：https://codeforces.com/blog/entry/48868
    template <typename F>
    size_t extreme(const F &dir) const {
        const auto &p = this->p;
        const auto check = [&](const size_t i) { 
            return dir(p[i]).toleft(p[this->nxt(i)] - p[i]) >= 0; 
        };
        const auto dir0 = dir(p[0]);
        const auto check0 = check(0);
        if (!check0 && check(p.size() - 1)) { return 0; }
        const auto cmp = [&](const Point &v) {
            const size_t vi = &v - p.data();
            if (vi == 0) { return 1; }
            const auto checkv = check(vi);
            const auto t = dir0.toleft(v - p[0]);
            if (vi == 1 && checkv == check0 && t == 0) { return 1; }
            return checkv ^ (checkv == check0 && t <= 0);
        };
        return partition_point(p.begin(), p.end(), cmp) - p.begin();
    }

    // 过凸多边形外一点求凸多边形的切线，返回切点下标
    // 复杂度 O(logn)
    // 必须保证点在多边形外
    pair<size_t, size_t> tangent(const point<T> &a) const {
        const size_t i = extreme([&](const point<T> &u) { return u - a; });
        const size_t j = extreme([&](const point<T> &u) { return a - u; });
        return {i, j};
    }

    // 求平行于给定直线的凸多边形的切线，返回切点下标
    // 复杂度 O(logn)
    pair<size_t, size_t> tangent(const line<T> &a) const {
        const size_t i = extreme([&](...) { return a.v; });
        const size_t j = extreme([&](...) { return -a.v; });
        return {i, j};
    }
};

using Convex = convex<_T>;

// 点集的凸包
// Andrew 算法，复杂度 O(nlogn)
Convex convexhull(vector<Point> p) {
    vector<Point> st;
    if (p.empty()) { return Convex{st}; }
    sort(p.begin(), p.end());
    const auto check = [](const vector<Point> &st, const Point &u) {
        const auto back1 = st.back(), back2 = *prev(st.end(), 2);
        return (back1 - back2).toleft(u - back1) <= 0;
    };
    for (const Point &u : p) {
        while (st.size() > 1 && check(st, u)) { st.pop_back(); }
        st.push_back(u);
    }
    size_t k = st.size();
    p.pop_back();
    reverse(p.begin(), p.end());
    for (const Point &u : p) {
        while (st.size() > k && check(st, u)) { st.pop_back(); }
        st.push_back(u);
    }
    st.pop_back();
    return Convex{st};
}

signed main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n;
    cin >> n;

    vector<Point> p(n);
    for (int i = 0; i < n; i ++ ) {
        cin >> p[i];
    }

    cout << "凸包的周长为：" << fixed << setprecision(2) << convexhull(p).circ() << "\n";

    return 0;
}