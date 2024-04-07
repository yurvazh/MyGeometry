//Created by Yury Vazhenin
#include <iostream>
#include <vector>
#include <initializer_list>
#include <cmath>

class Point;
class Line;
class Ellipse;
class Circle;
class Rectangle;
class Polygon;
class Square;
class Triangle;

const long double eps = 0.0001;

bool small(double d) {
    return (d < eps && (d > -eps));
}

class Point {
public:
    double x, y;

    Point() = default;

    Point (double x, double y);

    bool operator==(const Point& other) const;

    void rotatePoint(const Point& center, double angle);

    void reflect (const Point& center);

    Point operator+(const Point& other) const;

    Point operator-(const Point& other) const;

    Point operator*(double coefficient) const;

    void scale(const Point& center, const double& coefficient);
};


Point::Point(double x, double y) : x(x), y(y) {}

Point pointReflect(const Point& p, const Point& center) {
    return Point(2.0 * center.x - p.x, 2.0 * center.y - p.y);
}

void Point::scale (const Point& center, const double& coefficient) {
    *this = center + (*this - center) * coefficient;
}

void Point::reflect(const Point& center) {
    *this = pointReflect(*this, center);
}

void Point::rotatePoint(const Point& center, double angle) {
    double x_new = center.x + (x - center.x) * cos(angle) - (y - center.y) * sin(angle);
    double y_new = center.y + (y - center.y) * cos(angle) + (x - center.x) * sin(angle);
    *this = {x_new, y_new};
}

double dist (const Point& first_point, const Point& second_point) {
    double dx = first_point.x - second_point.x;
    double dy = first_point.y - second_point.y;
    return sqrt(dx * dx + dy * dy);
}

double crossProduct (const Point& a, const Point& b) {
    return (a.x * b.y - a.y * b.x);
}

double scalarProduct (const Point& a, const Point& b) {
    return a.x * b.x + a.y * b.y;
}

double crossProduct (const Point& a, const Point& b, const Point& c) {
    Point ab(b.x - a.x, b.y - a.y);
    Point ac(c.x - a.x, c.y - a.y);
    return crossProduct(ab, ac);
}

Point Point::operator+(const Point& other) const {
    return Point(x + other.x, y + other.y);
}

Point Point::operator-(const Point& other) const {
    return Point(x - other.x, y - other.y);
}

Point Point::operator*(double coefficient) const {
    return Point(x * coefficient, y * coefficient);
}

Point divide(const Point& a, const Point& b, double coefficient) {
    return a + (b - a) * (coefficient / (1.0 + coefficient));
}

bool Point::operator==(const Point& other) const {
    return small(x - other.x) && small(y - other.y);
}

Point pointScale(const Point& p, const Point& center, double coefficient) {
    Point answer = p;
    answer.scale(center, coefficient);
    return answer;
}

Point pointRotate (const Point& p, const Point& center, double angle) {
    Point answer = p;
    answer.rotatePoint(center, angle);
    return answer;
}

class Line {
public:
    double a, b, c;
public:
    Line (const Point& first_point, const Point& second_point);

    Line (double coefficient, double shift);

    Line (const Point& single_point, double coefficient);

    bool operator==(const Line& other) const;

    Point intersection(const Line& l) const;

    double distanceToPoint(const Point& point) const;
};

Line::Line(const Point& first_point, const Point& second_point) {
    a = first_point.y - second_point.y;
    b = second_point.x - first_point.x;
    c = - a * first_point.x - b * first_point.y;
}

Line::Line(double coefficient, double shift) {
    a = coefficient;
    b = -1.0;
    c = shift;
}

Line::Line(const Point& point, double coefficient) {
    a = -coefficient;
    b = 1.0;
    c = - (a * point.x + b * point.y);
}

bool Line::operator==(const Line& other) const {
    return (small(a * other.b - b * other.a)) && (small(c * other.b - b * other.c));
}

Point Line::intersection (const Line& other_line) const {
    Point answer;
    answer.x = (b * other_line.c - c * other_line.b) / (a * other_line.b - b * other_line.a);
    answer.y = (c * other_line.a - a * other_line.c) / (a * other_line.b - b * other_line.a);
    return answer;
}

double Line::distanceToPoint(const Point &point) const {
    Line orthogonal_line (point, point + Point(a, b));
    return dist(point, intersection(orthogonal_line));
}

Point pointReflect (const Point& p, const Line& axis) {
    Line orthogonal_line (p, Point(p.x + axis.a, p.y + axis.b));
    Point center = axis.intersection(orthogonal_line);
    return pointReflect(p, center);
}

class Shape {
public:
    virtual double perimeter() const = 0;

    virtual double area() const = 0;

    virtual bool containsPoint(const Point& point) const = 0;

    virtual void rotate(const Point& center, double angle) = 0;

    virtual void reflect(const Point& center) = 0;

    virtual void reflect(const Line& axis) = 0;

    virtual void scale(const Point& center, double coefficient) = 0;

    virtual bool operator==(const Shape& other) const = 0;

    virtual bool isCongruentTo(const Shape& other) const = 0;

    virtual bool isSimilarTo(const Shape& other) const = 0;

    virtual ~Shape() {};
};



class Polygon : public Shape {
protected:
    std::vector<Point> vertices;
public:
    Polygon() {};

    Polygon(const std::vector<Point>& vertices_vector);

    Polygon(std::initializer_list<Point> vertices_list);

    template <typename... Points>
    Polygon(Points... points) : vertices({points...}) {}

    int verticesCount() const;

    std::vector<Point> getVertices() const;

    bool isConvex() const;

    bool containsPoint(const Point& point) const override;

    double perimeter() const override;

    double area() const override;

    void scale(const Point& center, double coefficient) override;

    void reflect(const Point& center) override;

    void reflect(const Line& axis) override;

    void rotate (const Point& center, double angle) override;

    bool operator==(const Shape& other) const override;

    bool operator==(const Polygon& other) const;

    bool isSimilarTo(const Shape& other) const override;

    bool isCongruentTo(const Shape& other) const override;
};


Polygon::Polygon(const std::vector<Point>& vertices) : vertices(vertices) {}

Polygon::Polygon(std::initializer_list<Point> vertices) : vertices(vertices) {}

int Polygon::verticesCount() const {
    return vertices.size();
}

void Polygon::rotate(const Point& center, double angle) {
    for (int i = 0; i < verticesCount(); i++) {
        vertices[i].rotatePoint(center, angle);
    }
}

std::vector<Point> Polygon::getVertices() const {
    return vertices;
}

bool Polygon::containsPoint(const Point &point) const {
    double angle_sum = 0.0;
    for (int i = 0; i < verticesCount(); i++) {
        const Point& previous_point = vertices[i];
        const Point& next_point = vertices[(i + 1) % verticesCount()];
        double cross = crossProduct((previous_point - point), (next_point - point));
        double scalar = scalarProduct((previous_point - point), (next_point - point));
        if (point == previous_point)
            return true;
        angle_sum += atan2(cross, scalar);
    }
    return (angle_sum > M_PI || angle_sum < -M_PI);
}

double Polygon::perimeter() const {
    double result = 0.0;
    for (int i = 0; i < verticesCount(); i++)
        result += dist(vertices[i], vertices[(i + 1) % verticesCount()]);
    return result;
}

double Polygon::area() const {
    double result = 0.0;
    for (int i = 0; i < verticesCount(); i++) {
        Point first_point = vertices[i];
        Point second_point = vertices[(i + 1) % verticesCount()];
        result += (second_point.x - first_point.x) * (first_point.y + second_point.y);
    }
    result /= 2.0;
    return std::max(result, 0.0 - result);
}

bool Polygon::isConvex() const {
    int count_negative_angles = 0, count_positive_angles = 0;
    for (int i = 0; i < verticesCount(); i++) {
        double cross = crossProduct(vertices[i], vertices[(i + 1) % verticesCount()], vertices[(i + 2) % verticesCount()]);
        if (cross < 0.0)
            ++count_negative_angles;
        else
            ++count_positive_angles;
    }
    return (count_positive_angles == verticesCount()) || (count_negative_angles == verticesCount());
}

void Polygon::scale (const Point& center, double coefficient) {
    for (int i = 0; i < verticesCount(); i++) {
        vertices[i].scale(center, coefficient);
    }
}

void Polygon::reflect(const Point& center) {
    scale(center, -1.0);
}

void Polygon::reflect(const Line& axis) {
    for (int i = 0; i < verticesCount(); i++) {
        vertices[i] = pointReflect(vertices[i], axis);
    }
}

bool Polygon::operator==(const Shape& other) const {
    const Polygon* other_polygon = dynamic_cast<const Polygon*>(&other);
    if (other_polygon == nullptr)
        return false;
    if (verticesCount() != other_polygon->verticesCount())
        return false;
    std::vector<Point> firstVertices = vertices;
    const std::vector<Point>& secondVertices = other_polygon->vertices;
    int in = 0;
    for (int i = 0; i < verticesCount(); i++) {
        if (firstVertices[i] == secondVertices[0]) {
            in = i;
            break;
        }
    }
    std::rotate(firstVertices.begin(), firstVertices.begin() + in, firstVertices.end());
    bool equal = true;
    for (int i = 0; i < verticesCount(); i++) {
        equal &= (firstVertices[i] == secondVertices[i]);
        if (!equal) break;
    }
    if (equal)
        return true;
    std::reverse(firstVertices.begin() + 1, firstVertices.end());
    for (int i = 0; i < verticesCount(); i++) {
        if (firstVertices[i] != secondVertices[i]) {
            return false;
        }
    }
    return true;
}

bool check_similarity (const std::vector<Point>& a, const std::vector<Point>& b) {
    int n = a.size();
    int count_positive_angles = 0, count_negative_angles = 0;
    for (int i = 0; i < n; i++) {
        double cross_a_vertices = crossProduct(a[i] - a[(i + 1) % n], a[(i + 2) % n] - a[(i + 1) % n]);
        double scalar_a_vertices = scalarProduct(a[i] - a[(i + 1) % n], a[(i + 2) % n] - a[(i + 1) % n]);
        double first_polygon_angle = atan2(cross_a_vertices, scalar_a_vertices);
        double cross_b_vertices = crossProduct(b[i] - b[(i + 1) % n], b[(i + 2) % n] - b[(i + 1) % n]);
        double scalar_b_vertices = scalarProduct(b[i] - b[(i + 1) % n], b[(i + 2) % n] - b[(i + 1) % n]);
        double second_polygon_angle = atan2(cross_b_vertices, scalar_b_vertices);
        if (small(first_polygon_angle - second_polygon_angle))
            count_positive_angles++;
        else if (small(first_polygon_angle + second_polygon_angle))
            count_negative_angles++;
        else
            return false;
        double firstPolygonFraction = dist(a[(i + 1) % n], a[i]) / dist(b[(i + 1) % n], b[i]);
        double secondPolygonFraction = dist(a[(i + 2) % n], a[(i + 1) % n]) / dist(b[(i + 2) % n], b[(i + 1) % n]);
        if (!small(firstPolygonFraction - secondPolygonFraction))
            return false;
    }
    return (count_negative_angles == n || count_positive_angles == n);
}

bool Polygon::isSimilarTo(const Shape& other) const {
    const Polygon* other_polygon = dynamic_cast<const Polygon*>(&other);
    if (other_polygon == nullptr)
        return false;
    if (verticesCount() != other_polygon->verticesCount())
        return false;
    std::vector<Point> other_vertices = other_polygon->vertices;
    for (int i = 0; i < verticesCount(); i++) {
        if (check_similarity(vertices, other_vertices))
            return true;
        std::rotate(other_vertices.begin(), other_vertices.begin() + 1, other_vertices.end());
    }
    std::reverse(other_vertices.begin() + 1, other_vertices.end());
    for (int i = 0; i < verticesCount(); i++) {
        if (check_similarity(vertices, other_vertices))
            return true;
        std::rotate(other_vertices.begin(), other_vertices.begin() + 1, other_vertices.end());
    }
    return false;
}

bool Polygon::isCongruentTo(const Shape &other) const {
    const Polygon* other_polygon = dynamic_cast<const Polygon*>(&other);
    if (other_polygon == nullptr)
        return false;
    return (small(perimeter() - other_polygon->perimeter()) && isSimilarTo(other));
}

bool Polygon::operator==(const Polygon& other) const {
    return operator==(static_cast<const Shape&>(other));
}



class Ellipse : public Shape {
protected:
    Point f1, f2;
    double a;
public:
    Ellipse(const Point& f1, const Point& f2, double a);

    std::pair<Point, Point> focuses() const;

    Point center() const;

    double eccentricity() const;

    std::pair<Line, Line> directrices() const;

    double perimeter() const override;

    double area() const override;

    void reflect(const Point& center) override;

    void reflect (const Line& axis) override;

    void scale (const Point& center, double coefficient) override;

    void rotate(const Point& center, double angle) override;

    bool containsPoint(const Point& point) const override;

    bool operator==(const Shape& other) const override;

    bool isCongruentTo(const Shape& other) const override;

    bool isSimilarTo(const Shape& other) const override;
};


Ellipse::Ellipse(const Point& f1, const Point& f2, double a) : f1(f1), f2(f2), a(a / 2.0) {}

std::pair<Point, Point> Ellipse::focuses() const {
    return std::make_pair(f1, f2);
}

Point Ellipse::center() const {
    return Point((f1.x + f2.x) / 2.0, (f1.y + f2.y) / 2.0);
}

double Ellipse::eccentricity() const {
    return dist(f1, f2) / (2.0 * a);
}

std::pair<Line, Line> Ellipse::directrices() const {
    Point ellipseCenter = center();
    Point firstDirectrixPoint = pointScale(f2, ellipseCenter, (1.0 / (eccentricity() * eccentricity())));
    Point secondDirectrixPoint = pointScale(f1, ellipseCenter, (1.0 / (eccentricity() * eccentricity())));
    Line l (f1, f2);
    Line first_directrix(firstDirectrixPoint, Point(firstDirectrixPoint.x + l.a, firstDirectrixPoint.y + l.b));
    Line second_directrix(secondDirectrixPoint, Point(secondDirectrixPoint.x + l.a, secondDirectrixPoint.y + l.b));
    return std::make_pair(first_directrix, second_directrix);
}

void Ellipse::reflect(const Point& center) {
    f1.reflect(center);
    f2.reflect(center);
}

void Ellipse::reflect(const Line& axis) {
    f1 = pointReflect(f1, axis);
    f2 = pointReflect(f2, axis);
}

void Ellipse::scale(const Point& center, double coefficient) {
    f1.scale(center, coefficient);
    f2.scale(center, coefficient);
    a *= coefficient;
}

double Ellipse::perimeter() const {
    double b = sqrt(a * a * (1.0 - eccentricity() * eccentricity()));
    return M_PI * (3.0 * (a + b) - sqrt((3.0 * a + b) * (a + 3.0 * b)));
}

double Ellipse::area() const {
    double b = sqrt(a * a * (1.0 - eccentricity() * eccentricity()));
    return M_PI * a * b;
}

void Ellipse::rotate(const Point& center, double angle) {
    f1.rotatePoint(center, angle);
    f2.rotatePoint(center, angle);
}

bool Ellipse::containsPoint(const Point &point) const {
    double distanceDifference = (dist(point, f1) + dist(point, f2) - 2.0 * a);
    return (distanceDifference <= 0.0);
}

bool Ellipse::operator==(const Shape& other) const {
    const Ellipse* other_ellipse = dynamic_cast<const Ellipse*>(&other);
    if (other_ellipse == nullptr)
        return false;
    std::pair<Point, Point> otherFocuses = other_ellipse->focuses();
    if (otherFocuses == focuses() || (otherFocuses == std::make_pair(f2, f1)))
        return small(a - other_ellipse->a);
    return false;
}

bool Ellipse::isCongruentTo(const Shape& other) const {
    const Ellipse* other_ellipse = dynamic_cast<const Ellipse*>(&other);
    if (other_ellipse == nullptr)
        return false;
    std::pair<Point, Point> otherFocuses = other_ellipse->focuses();
    double distanceBetweenFirstFocuses = dist(f1, f2);
    double distanceBetweenSecondFocuses = dist(otherFocuses.first, otherFocuses.second);
    if (small(distanceBetweenFirstFocuses - distanceBetweenSecondFocuses))
        return a == other_ellipse->a;
    return false;
}

bool Ellipse::isSimilarTo(const Shape& other) const {
    const Ellipse* other_ellipse = dynamic_cast<const Ellipse*>(&other);
    if (other_ellipse == nullptr)
        return false;
    return eccentricity() == other_ellipse->eccentricity();
}

class Circle : public Ellipse {
public:
    Circle (const Point& center, double radius);

    double radius() const;

    Point center() const;

    bool operator==(const Circle& other) const;
};

Circle::Circle(const Point& center, double radius) : Ellipse(center, center, 2.0 * radius) {}

double Circle::radius() const {
    return a;
}

Point Circle::center() const {
    return f1;
}

bool Circle::operator==(const Circle& other) const {
    return (center() == other.center()) && (radius() == other.radius());
}


class Rectangle : public Polygon {
public:
    Rectangle (const Point& p1, const Point& p2, double ratio);

    Point center() const;

    std::pair<Line, Line> diagonalis() const;
};


Rectangle::Rectangle(const Point& firstVertex, const Point& secondVertex, double ratio) {
    if (ratio < 1.0) {
        ratio = 1.0 / ratio;
    }
    double cos_a = 1.0 / (1.0 + ratio * ratio);
    double sin_a = 1.0 / (1.0 + (1.0 / (ratio * ratio)));
    Line firstSide {firstVertex, pointRotate(secondVertex, firstVertex, atan2(sin_a, cos_a))};
    std::swap(sin_a, cos_a);
    sin_a = - sin_a;
    Line secondSide {secondVertex, pointRotate(firstVertex, secondVertex, atan2(sin_a, cos_a))};
    Point thirdVertex = firstSide.intersection(secondSide);
    Point fourthVertex = pointReflect(thirdVertex, Point((firstVertex.x + secondVertex.x) / 2.0, (firstVertex.y + secondVertex.y) / 2.0));
    vertices = {firstVertex, thirdVertex, secondVertex, fourthVertex};
}

Point Rectangle::center() const {
    return Point((vertices[0].x + vertices[2].x) / 2.0, (vertices[0].y + vertices[2].y) / 2.0);
}

std::pair<Line, Line> Rectangle::diagonalis() const {
    return std::make_pair(Line(vertices[0], vertices[2]), Line(vertices[1], vertices[3]));
}

class Square : public Rectangle {
    double a;
public:
    Square (const Point& firstVertex, const Point& secondVertex);

    virtual double perimeter() const override;

    virtual double area() const override;

    Circle circumscribedCircle() const;

    Circle inscribedCircle() const;
};


Square::Square (const Point& firstVertex, const Point& secondVertex) : Rectangle(firstVertex, secondVertex, 1.0)
    , a(dist(firstVertex, secondVertex) / sqrt(2.0)) {}

double Square::perimeter() const {
    return 4.0 * a;
}

double Square::area() const {
    return a * a;
}

Circle Square::circumscribedCircle() const {
    return Circle(center(), a / 2.0);
}

Circle Square::inscribedCircle() const {
    return Circle(center(), a / sqrt(2.0));
}

class Triangle : public Polygon {
public:
    using Polygon::Polygon;

    Circle circumscribedCircle() const;

    Circle inscribedCircle() const;

    Point centroid() const;

    Point orthocenter() const;

    Line EulerLine() const;

    Circle ninePointsCircle() const;
};

Point Triangle::centroid() const {
    return (vertices[0] + vertices[1] + vertices[2]) * (1.0 / 3.0);
}

Point Triangle::orthocenter() const {
    const Point& a = vertices[0];
    const Point& b = vertices[1];
    const Point& c = vertices[2];
    Line bc(b, c);
    Line ac(a, c);
    Line altitudeA(a, a + Point(bc.a, bc.b));
    Line altitudeB(b, b + Point(ac.a, ac.b));
    return altitudeA.intersection(altitudeB);
}

Circle Triangle::circumscribedCircle() const {
    const Point& a = vertices[0];
    const Point& b = vertices[1];
    const Point& c = vertices[2];
    Line bc (b, c);
    Line ab(a, b);
    Point bcCenter = (b + c) * 0.5;
    Point abCenter = (a + b) * 0.5;
    Line bcPerpendicular (bcCenter, bcCenter + Point(bc.a, bc.b));
    Line abPerpendicular (abCenter, abCenter + Point(ab.a, ab.b));
    Point circleCenter = abPerpendicular.intersection(bcPerpendicular);
    return Circle(circleCenter, dist(circleCenter, a));
}

Circle Triangle::ninePointsCircle() const {
    const Point& abCenter = (vertices[0] + vertices[1]) * 0.5;
    const Point& bcCenter = (vertices[0] + vertices[2]) * 0.5;
    const Point& acCenter = (vertices[2] + vertices[1]) * 0.5;
    return Triangle({abCenter, bcCenter, acCenter}).circumscribedCircle();
}

Line Triangle::EulerLine() const {
    return Line(centroid(), orthocenter());
}

Circle Triangle::inscribedCircle() const {
    const Point& a = vertices[0];
    const Point& b = vertices[1];
    const Point& c = vertices[2];
    Line bisectrixA (a, divide(b, c, dist(a, b) / dist(a, c)));
    Line bisectrixB (b, divide(a, c, dist(a, b) / dist(b, c)));
    Point circleCenter = bisectrixA.intersection(bisectrixB);
    return Circle(circleCenter, Line(a, b).distanceToPoint(circleCenter));
}
