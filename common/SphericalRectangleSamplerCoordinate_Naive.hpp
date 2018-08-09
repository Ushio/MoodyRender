	struct Rectangle {
		Rectangle() {}
		Rectangle(glm::dvec3 s, glm::dvec3 ex, glm::dvec3 ey)
		:_s(s)
		,_ex(ex)
		,_ey(ey) {
			_exLength = glm::length(ex);
			_eyLength = glm::length(ey);
			_x = ex / _exLength;
			_y = ey / _exLength;
			_z = glm::cross(_x, _y);
		}
		glm::dvec3 sample(double u, double v) const {
			return _s + _ex * u + _ey * v;
		}
		double area() const {
			return _exLength * _eyLength;
		}
		glm::dvec3 normal() const {
			return _z;
		}
		glm::dvec3 s() const {
			return _s;
		}
		glm::dvec3 ex() const {
			return _ex;
		}
		glm::dvec3 ey() const {
			return _ey;
		}
		glm::dvec3 x() const {
			return _x;
		}
		glm::dvec3 y() const {
			return _y;
		}
		glm::dvec3 z() const {
			return _z;
		}
		double exLength() const {
			return _exLength;
		}
		double eyLength() const {
			return _eyLength;
		}
		glm::dvec3 _s;
		glm::dvec3 _ex;
		glm::dvec3 _ey;
		double _exLength = 0.0;
		double _eyLength = 0.0;
		glm::dvec3 _x;
		glm::dvec3 _y;
		glm::dvec3 _z;
	};

	class SphericalRectangleSamplerCoordinate {
	public:
		SphericalRectangleSamplerCoordinate(const Rectangle &rectangle, const glm::dvec3 &o):_rectangle(rectangle) {
			_o = o;

			glm::dvec3 d = rectangle.s() - o;
			_x0 = glm::dot(d, rectangle.x());
			_y0 = glm::dot(d, rectangle.y());
			_z0 = glm::dot(d, rectangle.z());

			_x1 = _x0 + rectangle.exLength();
			_y1 = _y0 + rectangle.eyLength();

			// z flip
			if (_z0 > 0.0) {
				_z0 *= -1.0;
				_rectangle._z *= -1.0;
			}

			_v[0][0] = glm::vec3(_x0, _y0, _z0);
			_v[0][1] = glm::vec3(_x0, _y1, _z0);
			_v[1][0] = glm::vec3(_x1, _y0, _z0);
			_v[1][1] = glm::vec3(_x1, _y1, _z0);

			auto calculate_n = [](glm::vec3 a, glm::vec3 b) {
				glm::vec3 c = glm::cross(a, b);
				return glm::normalize(c);
			};

			_n[0] = calculate_n(_v[0][0], _v[1][0]);
			_n[1] = calculate_n(_v[1][0], _v[1][1]);
			_n[2] = calculate_n(_v[1][1], _v[0][1]);
			_n[3] = calculate_n(_v[0][1], _v[0][0]);

			_g[0] = std::acos(glm::dot(-_n[0], _n[1]));
			_g[1] = std::acos(glm::dot(-_n[1], _n[2]));
			_g[2] = std::acos(glm::dot(-_n[2], _n[3]));
			_g[3] = std::acos(glm::dot(-_n[3], _n[0]));
		}

		double solidAngle() const {
			return _g[0] + _g[1] + _g[2] + _g[3] - glm::two_pi<double>();
		}

		double minPhiU() const {
			return -_g[2] - _g[3] + glm::two_pi<double>();
		}
		double maxPhiU() const {
			return solidAngle() - _g[2] - _g[3] + glm::two_pi<double>();
		}
		glm::dvec3 sample(double u, double v) const {
			double AQ = solidAngle();

			double phi_u = u * AQ - _g[2] - _g[3] + glm::two_pi<double>();

			auto safeSqrt = [](double x) {
				return std::sqrt(std::max(x, 0.0));
			};

			double b0 = _n[0].z;
			double b1 = _n[2].z;
			double fu = (std::cos(phi_u) * b0 - b1) / std::sin(phi_u);
			double cu = std::copysign(1.0, fu) / safeSqrt(fu * fu + b0 * b0);
			double xu = -cu * _z0 / safeSqrt(1.0 - cu * cu);

			double d = std::sqrt(xu * xu + _z0 * _z0);
			double h0 = _y0 / std::sqrt(d * d + _y0 * _y0);
			double h1 = _y1 / std::sqrt(d * d + _y1 * _y1);
			double hv = glm::mix(h0, h1, v);
			double yv = hv * d / safeSqrt(1.0 - hv * hv);
			return _o + xu * _rectangle.x() + yv * _rectangle.y() + _z0 * _rectangle.z();
		}

		Rectangle _rectangle;

		glm::dvec3 _o;

		// local coord : 
		double _x0;
		double _x1;
		double _y0;
		double _y1;
		double _z0;

		glm::dvec3 _v[2][2];
		glm::dvec3 _n[4];
		double _g[4];
	};