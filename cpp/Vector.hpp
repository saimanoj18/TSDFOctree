#pragma once
#include <stdio.h>
#include <Eigen/Dense>
namespace BPVector{
    struct Vector3{
        Vector3(){
        }
        Vector3(double x, double y, double z){
            X = x;
            Y = y;
            Z = z;
        }

        void operator+=(const Vector3 &rhs){
            X += rhs.X;
            Y += rhs.Y;
            Z += rhs.Z;
        }

        void operator-=(const Vector3 &rhs){
            X -= rhs.X;
            Y -= rhs.Y;
            Z -= rhs.Z;
        }
        Vector3 operator-(const Vector3 &rhs){
            return Vector3( X - rhs.X, Y - rhs.Y, Z - rhs.Z);
        }
        Vector3 operator-(){
            return Vector3( -X, -Y, -Z);
        }
        Vector3 operator+(const Vector3 &rhs){
            return Vector3( X + rhs.X, Y + rhs.Y, Z + rhs.Z);
        }
        Vector3 operator*(const double scale){
            return Vector3(X * scale, Y * scale, Z * scale);
        }
        Vector3 operator^(const Vector3 &rhs){
            return Vector3(
                this->Y * rhs.Z - this->Z * rhs.Y,
                this->Z * rhs.X - this->X * rhs.Z,
                this->X * rhs.Y - this->Y * rhs.X
                );
        }
        double dot(const Vector3 &rhs){
            return X*rhs.X + Y*rhs.Y + Z*rhs.Z;
        }

        double norm(){
            return sqrt( X*X + Y*Y + Z*Z );
        }
        void normalize(){
            double Norm = norm();
            X /= Norm;
            Y /= Norm;
            Z /= Norm;
        }
        void print(){
            printf("<%lf %lf %lf>\t", X, Y, Z);
        }

        double X, Y, Z;
    };

    struct Vector4{
        Vector4(){
            W = 0;
            X = 0;
            Y = 0;
            Z = 0;
        }
        Vector4(double w, double x, double y, double z){
            W = w;
            X = x;
            Y = y;
            Z = z;
        }

        void operator+=(const Vector4 &rhs){
            W += rhs.W;
            X += rhs.X;
            Y += rhs.Y;
            Z += rhs.Z;
        }

        void operator-=(const Vector4 &rhs){
            W -= rhs.W;
            X -= rhs.X;
            Y -= rhs.Y;
            Z -= rhs.Z;
        }
        Vector4 operator-(const Vector4 &rhs){
            return Vector4( W - rhs.W, X - rhs.X, Y - rhs.Y, Z - rhs.Z);
        }
        Vector4 operator+(const Vector4 &rhs){
            return Vector4( W + rhs.W, X + rhs.X, Y + rhs.Y, Z + rhs.Z);
        }
        Vector4 operator*(const double scale){
            return Vector4( W * scale, X * scale, Y * scale, Z * scale);
        }
        double norm() const {
            return sqrt(W*W + Y*Y + X*X + Z*Z);
        }
        void print(){
            printf("<%lf %lf %lf %lf>\n", W, X, Y, Z);
        }

        double W, X, Y, Z;
    };

    struct Quaternion : Vector4{
        Quaternion(){
            W = 1;
            X = 0;
            Y = 0;
            Z = 0;
        }
        Quaternion(double w, double x, double y, double z){
            W = w;
            X = x;
            Y = y;
            Z = z;
        }
        Quaternion(double angle, Vector3 axis){
            W = cos(angle/2);
            X = axis.X * sin(angle/2);
            Y = axis.Y * sin(angle/2);
            Z = axis.Z * sin(angle/2);
        }
        Quaternion(Vector3 Euler){
            W = cos(Euler.X/2)*cos(Euler.Y/2)*cos(Euler.Z/2) + sin(Euler.X/2)*sin(Euler.Y/2)*sin(Euler.Z/2);
            X = sin(Euler.X/2)*cos(Euler.Y/2)*cos(Euler.Z/2) - cos(Euler.X/2)*sin(Euler.Y/2)*sin(Euler.Z/2);
            Y = cos(Euler.X/2)*sin(Euler.Y/2)*cos(Euler.Z/2) + sin(Euler.X/2)*cos(Euler.Y/2)*sin(Euler.Z/2);
            Z = cos(Euler.X/2)*cos(Euler.Y/2)*sin(Euler.Z/2) - sin(Euler.X/2)*sin(Euler.Y/2)*cos(Euler.Z/2);
        }
        Quaternion(Eigen::Matrix<float, 3, 3, Eigen::RowMajor> Rotation){
            float tr = Rotation(0,0) + Rotation(1,1) + Rotation(2,2);

            if (tr > 0) {
              float S = sqrt(tr+1.0) * 2; // S=4*W
              W = 0.25 * S;
              X = (Rotation(2,1) - Rotation(1,2)) / S;
              Y = (Rotation(0,2) - Rotation(2,0)) / S;
              Z = (Rotation(1,0) - Rotation(0,1)) / S;
            } else if ((Rotation(0,0) > Rotation(1,1))&(Rotation(0,0) > Rotation(2,2))) {
              float S = sqrt(1.0 + Rotation(0,0) - Rotation(1,1) - Rotation(2,2)) * 2; // S=4*X
              W = (Rotation(2,1) - Rotation(1,2)) / S;
              X = 0.25 * S;
              Y = (Rotation(0,1) + Rotation(1,0)) / S;
              Z = (Rotation(0,2) + Rotation(2,0)) / S;
            } else if (Rotation(1,1) > Rotation(2,2)) {
              float S = sqrt(1.0 + Rotation(1,1) - Rotation(0,0) - Rotation(2,2)) * 2; // S=4*Y
              W = (Rotation(0,2) - Rotation(2,0)) / S;
              X = (Rotation(0,1) + Rotation(1,0)) / S;
              Y = 0.25 * S;
              Z = (Rotation(1,2) + Rotation(2,1)) / S;
            } else {
              float S = sqrt(1.0 + Rotation(2,2) - Rotation(0,0) - Rotation(1,1)) * 2; // S=4*Z
              W = (Rotation(1,0) - Rotation(0,1)) / S;
              X = (Rotation(0,2) + Rotation(2,0)) / S;
              Y = (Rotation(1,2) + Rotation(2,1)) / S;
              Z = 0.25 * S;
            }

        }

        Quaternion operator-(const Quaternion &rhs){
            return Quaternion( W - rhs.W, X - rhs.X, Y - rhs.Y, Z - rhs.Z);
        }
        Quaternion operator+(const Quaternion &rhs){
            return Quaternion( W + rhs.W, X + rhs.X, Y + rhs.Y, Z + rhs.Z);
        }
        Quaternion operator*(const double scale){
            return Quaternion( W * scale, X * scale, Y * scale, Z * scale);
        }
        Quaternion operator*(const Quaternion &rhs){
            return Quaternion(
                        W * rhs.W - X * rhs.X - Y * rhs.Y - Z * rhs.Z,
                        W * rhs.X + X * rhs.W + Y * rhs.Z - Z * rhs.Y,
                        W * rhs.Y - X * rhs.Z + Y * rhs.W + Z * rhs.X,
                        W * rhs.Z + X * rhs.Y - Y * rhs.X + Z * rhs.W
                    );
        }
        Vector3 EulerAngles(Vector3 direction = Vector3(0,0,0) ) const {
            Vector3 angles(
                        atan2( 2*(W*X + Y*Z), 1 - 2*(X*X + Y*Y) ),
                        asin( 2*(W*Y - Z*X) ),
                        atan2( 2*(W*Z + X*Y), 1 - 2*(Y*Y + Z*Z) )
                    );

            if(angles.norm() != 0){
                angles.X += direction.X * angles.X < 0 ? direction.X * 2.0*3.14159 : 0;
                angles.Y += direction.Y * angles.Y < 0 ? direction.Y * 2.0*3.14159 : 0;
                angles.Z += direction.Z * angles.Z < 0 ? direction.Z * 2.0*3.14159 : 0;
            }else{
                angles.X += direction.X * angles.X < 0 ? direction.X * 2.0*3.14159 : 0;
                angles.Y += direction.Y * angles.Y < 0 ? direction.Y * 2.0*3.14159 : 0;
                angles.Z += direction.Z * angles.Z < 0 ? direction.Z * 2.0*3.14159 : 0;
            }

            return angles;
        }

        void normalize(){
            double magnitude = norm();
            *this = this->operator*(1.0/magnitude);
        }
        Quaternion inverse() const {
            double div = norm();
            div *= div;

            return Quaternion(W/div,-X/div,-Y/div,-Z/div);
        }

        Vector3 Rotate(const Vector3 &v) const {
            Quaternion p(0, v.X, v.Y, v.Z);
            Quaternion ans = this->inverse() * p * (*this);

            return Vector3( ans.X, ans.Y, ans.Z);
            /*return Vector3(
                v.X * (1 - 2*Y*Y - 2*Z*Z) + v.Y * (2*X*Y - 2*W*Z) + v.Z * (2*X*Z + 2*W*Y),
                v.X * (2*X*Y + 2*W*Z) + v.Y * (1 - 2*X*X - 2*Z*Z) + v.Z * (2*Y*Z - 2*W*X),
                v.X * (2*X*Z - 2*W*Y) + v.Y * (2*Y*Z + 2*W*X) + v.Z * (1 - 2*X*X - 2*Y*Y)
                );*/
        }
        Eigen::Matrix<float, 3, 3, Eigen::RowMajor> RotationMatrix() const {
            Eigen::Matrix<float, 3, 3, Eigen::RowMajor> R;

            R(0,0) = (1 - 2*Y*Y - 2*Z*Z);	R(0,1) = (2*X*Y - 2*W*Z);		R(0,2) = (2*X*Z + 2*W*Y);
            R(1,0) = (2*X*Y + 2*W*Z);		R(1,1) =  (1 - 2*X*X - 2*Z*Z);	R(1,2) = (2*Y*Z - 2*W*X);
            R(2,0) = (2*X*Z - 2*W*Y);		R(2,1) = (2*Y*Z + 2*W*X);		R(2,2) = (1 - 2*X*X - 2*Y*Y);

            return R;
        }


    };

    struct State3D{
        State3D(){
            this->position = Vector3(0,0,0);
            this->angle = Vector3(0,0,0);
        }
        State3D(Vector3 position, Vector3 angle){
            this->position = position;
            this->angle = angle;
        }

        void operator+=(const State3D &rhs){
            position += rhs.position;
            angle += rhs.angle;
        }

        void operator-=(const State3D &rhs){
            position -= rhs.position;
            angle -= rhs.angle;
        }
        State3D operator-(const State3D &rhs){
            return State3D(position - rhs.position, angle - rhs.angle);
        }
        State3D operator+(const State3D &rhs){
            return State3D( position + rhs.position, angle + rhs.angle);
        }
        State3D operator*(const double scale){
            return State3D( position * scale, angle * scale );
        }

        Vector3 position;
        Vector3 angle;
    };

    struct Motion3D{
        Motion3D(){
            position = Vector3(0,0,0);
            q = Quaternion(0,0,0,0);
        }
        Motion3D(Vector3 position, Quaternion q){
            this->position = position;
            this->q = q;
        }

        void operator+=(const Motion3D &rhs){
            position += rhs.position;
            q += rhs.q;
        }

        void operator-=(const Motion3D &rhs){
            position -= rhs.position;
            q -= rhs.q;
        }
        Motion3D operator-(const Motion3D &rhs){
            return Motion3D(position - rhs.position, q - rhs.q);
        }
        Motion3D operator+(const Motion3D &rhs){
            return Motion3D( position + rhs.position, q + rhs.q);
        }
        Motion3D operator*(const double scale){
            return Motion3D( position * scale, q * scale );
        }
        double norm(){
            return sqrt( position.norm() + q.norm() );
        }

        Vector3 position;
        Quaternion q;
    };
}
