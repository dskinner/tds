// Package tds provides methods for delaunay triangulation.
package tds

import (
	"fmt"
	"log"
	"math"
)

func sqrt(x float32) float32 {
	return float32(math.Sqrt(float64(x)))
}

type Vec2 [2]float32

// Sub stores the result of a-b in v.
func (v *Vec2) Sub(a, b *Vec2) {
	v[0] = a[0] - b[0]
	v[1] = a[1] - b[1]
}

// Add stores the result of a+b in v.
func (v *Vec2) Add(a, b *Vec2) {
	v[0] = a[0] + b[0]
	v[1] = a[1] + b[1]
}

func (v *Vec2) Dot(a *Vec2) float32 {
	return v[0]*a[0] + v[1]*a[1]
}

func (v *Vec2) DivScalar(a float32) {
	v[0] /= a
	v[1] /= a
}

type Vec3 [3]float32

// Sub stores the result of a-b in v.
func (v *Vec3) Sub(a, b *Vec3) {
	v[0] = a[0] - b[0]
	v[1] = a[1] - b[1]
	v[2] = a[2] - b[2]
}

// Add stores the result of a+b in v.
func (v *Vec3) Add(a, b *Vec3) {
	v[0] = a[0] + b[0]
	v[1] = a[1] + b[1]
	v[2] = a[2] + b[2]
}

func (v *Vec3) Dot(a *Vec3) float32 {
	return v[0]*a[0] + v[1]*a[1] + v[2]*a[2]
}

type Vec4 [4]float32

// Sub stores the result of a-b in v.
func (v *Vec4) Sub(a, b *Vec4) {
	v[0] = a[0] - b[0]
	v[1] = a[1] - b[1]
	v[2] = a[2] - b[2]
	v[3] = a[3] - b[3]
}

// Add stores the result of a+b in v.
func (v *Vec4) Add(a, b *Vec4) {
	v[0] = a[0] + b[0]
	v[1] = a[1] + b[1]
	v[2] = a[2] + b[2]
	v[3] = a[3] + b[3]
}

func (v *Vec4) Dot(a *Vec4) float32 {
	return v[0]*a[0] + v[1]*a[1] + v[2]*a[2] + v[3]*a[3]
}

type Mat2 [2]Vec2

// Det returns determinant of m.
func (m Mat2) Det() float32 {
	return m[0][0]*m[1][1] - m[0][1]*m[1][0]
}

type Mat3 [3]Vec3

// Det returns determinant of m.
func (m Mat3) Det() float32 {
	// [0,0  0,1  0,2]
	// [1,0  1,1  1,2]
	// [2,0  2,1  2,2]
	return m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1]) -
		m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0]) +
		m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0])
}

type Mat4 [4]Vec4

// Det returns determinant of m.
func (m Mat4) Det() float32 {
	// [0,0  0,1  0,2  0,3]
	// [1,0  1,1  1,2  1,3]
	// [2,0  2,1  2,2  2,3]
	// [3,0  3,1  3,2  3,3]
	return m[0][0]*(m[1][1]*(m[2][2]*m[3][3]-m[2][3]*m[3][2])-
		m[1][2]*(m[2][1]*m[3][3]-m[2][3]*m[3][1])+
		m[1][3]*(m[2][1]*m[3][2]-m[2][2]*m[3][1])) -
		m[0][1]*(m[1][0]*(m[2][2]*m[3][3]-m[2][3]*m[3][2])-
			m[1][2]*(m[2][0]*m[3][3]-m[2][3]*m[3][0])+
			m[1][3]*(m[2][0]*m[3][2]-m[2][2]*m[3][0])) +
		m[0][2]*(m[1][0]*(m[2][1]*m[3][3]-m[2][3]*m[3][1])-
			m[1][1]*(m[2][0]*m[3][3]-m[2][3]*m[3][0])+
			m[1][3]*(m[2][0]*m[3][1]-m[2][1]*m[3][0])) -
		m[0][3]*(m[1][0]*(m[2][1]*m[3][2]-m[2][2]*m[3][1])-
			m[1][1]*(m[2][0]*m[3][2]-m[2][2]*m[3][0])+
			m[1][2]*(m[2][0]*m[3][1]-m[2][1]*m[3][0]))
}

// Orient2D determines how points a, b, and c are arranged by the determinant
// of vectors a-c and b-c, returning a positive value if counter-clockwise,
// negative if clockwise, or zero if collinear.
func Orient2D(a, b, c Vec2) float32 {
	var m Mat2
	m[0].Sub(&a, &c)
	m[1].Sub(&b, &c)
	return m.Det()
}

// Orient3D determines how points a, b, c, and d are arranged by the determinant
// of vectors a-d, b-d, and c-d, returning a positive value if counter-clockwise,
// negative if clockwise, or zero if coplanar.
func Orient3D(a, b, c, d Vec3) float32 {
	var m Mat3
	m[0].Sub(&a, &d)
	m[1].Sub(&b, &d)
	m[2].Sub(&c, &d)
	return m.Det()
}

// InCircle returns a positive value if d lies inside the circle through
// a, b, and c oriented counter-clockwise, negative if outside, or zero if
// coplanar.
func InCircle(a, b, c, d *Vec2) float32 {
	// var v0, v1, v2 Vec2
	// v0.Sub(a, d)
	// v1.Sub(b, d)
	// v2.Sub(c, d)
	// return Mat3{
	// {v0[0], v0[1], v0.Dot(&v0)},
	// {v1[0], v1[1], v1.Dot(&v1)},
	// {v2[0], v2[1], v2.Dot(&v2)},
	// }.Det()
	m00 := a[0] - d[0]
	m01 := a[1] - d[1]
	m02 := m00*m00 + m01*m01
	m10 := b[0] - d[0]
	m11 := b[1] - d[1]
	m12 := m10*m10 + m11*m11
	m20 := c[0] - d[0]
	m21 := c[1] - d[1]
	m22 := m20*m20 + m21*m21
	return m00*(m11*m22-m12*m21) -
		m01*(m10*m22-m12*m20) +
		m02*(m10*m21-m11*m20)
}

// InSphere returns a positive value if d lies inside the sphere through
// a, b, c, and d oriented counter-clockwise, negative if outside, or zero if
// cospherical.
func InSphere(a, b, c, d, e Vec3) float32 {
	var v0, v1, v2, v3 Vec3
	v0.Sub(&a, &e)
	v1.Sub(&b, &e)
	v2.Sub(&c, &e)
	v3.Sub(&d, &e)
	return Mat4{
		{v0[0], v0[1], v0[2], v0.Dot(&v0)},
		{v1[0], v1[1], v1[2], v1.Dot(&v1)},
		{v2[0], v2[1], v2[2], v2.Dot(&v2)},
		{v3[0], v3[1], v3[2], v3.Dot(&v3)},
	}.Det()
}

type Store2D struct {
	m map[Mat2]Mat2
	g map[[3]Vec2][]Vec2
	v map[[3]Vec2]bool
}

func (st *Store2D) M() map[Mat2]Mat2 { return st.m }

func NewStore2D() *Store2D {
	return &Store2D{
		m: make(map[Mat2]Mat2),
		g: make(map[[3]Vec2][]Vec2),
		v: make(map[[3]Vec2]bool),
	}
}

func (st *Store2D) Vertices() (a []float32) {
	for mt := range st.m {
		a = append(a, mt[0][0], mt[0][1], mt[1][0], mt[1][1])
	}
	return
}

func (st *Store2D) AddTriangle(u, v, w Vec2) error {
	e0 := Mat2{u, v}
	if _, ok := st.m[e0]; ok {
		return fmt.Errorf("already contains edge %v", e0)
	}
	e1 := Mat2{v, w}
	if _, ok := st.m[e1]; ok {
		return fmt.Errorf("already contains edge %v", e1)
	}
	e2 := Mat2{w, u}
	if _, ok := st.m[e2]; ok {
		return fmt.Errorf("already contains edge %v", e2)
	}
	st.m[e0] = e1
	st.m[e1] = e2
	st.m[e2] = e0
	return nil
}

func (st *Store2D) DeleteTriangle(u, v, w Vec2) error {
	e0 := Mat2{u, v}
	if _, ok := st.m[e0]; !ok {
		return fmt.Errorf("does not contain edge %v", e0)
	}
	e1 := Mat2{v, w}
	if _, ok := st.m[e1]; !ok {
		return fmt.Errorf("does not contain edge %v", e1)
	}
	e2 := Mat2{w, u}
	if _, ok := st.m[e2]; !ok {
		return fmt.Errorf("does not contain edge %v", e2)
	}
	delete(st.m, e0)
	delete(st.m, e1)
	delete(st.m, e2)
	return nil
}

func (st *Store2D) Adjacent(u, v Vec2) (Vec2, bool) {
	m, ok := st.m[Mat2{u, v}]
	return m[1], ok
}

const weight = 0.001

// PointLocation returns one triangle whose open circumdisk contains u.
func (st *Store2D) PointLocation(u Vec2) (Vec2, Vec2, Vec2, error) {
	for e0, e1 := range st.m {
		v, w, x := e0[0], e1[0], e1[1]
		v3 := Vec3{v[0], v[1], v.Dot(&v) - weight}
		w3 := Vec3{w[0], w[1], w.Dot(&w) - weight}
		x3 := Vec3{x[0], x[1], x.Dot(&x) - weight}
		u3 := Vec3{u[0], u[1], u.Dot(&u) - weight}
		// log.Println(Orient3D(v3, w3, x3, u3))
		if Orient3D(v3, w3, x3, u3) > 0 {
			// if InCircle(&v, &w, &x, &u) > 0 {
			return v, w, x, nil
		}
	}
	return Vec2{}, Vec2{}, Vec2{}, fmt.Errorf("no triangle's open circumdisk contains %+v", u)
}

// var GhostVertex = Vec2{math.MaxFloat32, math.MaxFloat32}

var GhostVertex = Vec2{9, 9}

func (st *Store2D) InsertGhost() {
	var p [][3]Vec2
	for e0, e1 := range st.m {
		u, v, w := e0[0], e1[0], e1[1]
		if u == GhostVertex || v == GhostVertex || w == GhostVertex {
			continue
		}
		if _, ok := st.Adjacent(v, u); !ok {
			p = append(p, [3]Vec2{v, u, GhostVertex})
		}
		if _, ok := st.Adjacent(w, v); !ok {
			p = append(p, [3]Vec2{w, v, GhostVertex})
		}
		if _, ok := st.Adjacent(u, w); !ok {
			p = append(p, [3]Vec2{u, w, GhostVertex})
		}
	}
	for _, uvw := range p {
		st.AddTriangle(uvw[0], uvw[1], uvw[2])
	}
}

func (st *Store2D) dig(u Vec2, v, w Vec2) error {
	x, ok := st.Adjacent(w, v) // w, v, x is triangle on other side of edge v, w from u
	if !ok {
		return nil // triangle already deleted
	}
	u3 := Vec3{u[0], u[1], u.Dot(&u) - weight}
	v3 := Vec3{v[0], v[1], v.Dot(&v) - weight}
	w3 := Vec3{w[0], w[1], w.Dot(&w) - weight}
	x3 := Vec3{x[0], x[1], x.Dot(&x) - weight}
	if Orient3D(u3, v3, w3, x3) > 0 {
		// if InCircle(&u, &v, &w, &x) > 0 {
		if err := st.DeleteTriangle(w, v, x); err != nil {
			return err
		}
		if err := st.dig(u, v, x); err != nil {
			return err
		}
		if err := st.dig(u, x, w); err != nil {
			return err
		}
	} else {
		if err := st.AddTriangle(u, v, w); err != nil {
			return err
		}
	}
	return nil
}

func (st *Store2D) InsertVertex(u Vec2, v, w, x Vec2) error {
	if err := st.DeleteTriangle(v, w, x); err != nil {
		return err
	}
	if err := st.dig(u, v, w); err != nil {
		return err
	}
	if err := st.dig(u, w, x); err != nil {
		return err
	}
	if err := st.dig(u, x, v); err != nil {
		return err
	}
	st.InsertGhost()
	// st.FixAngles()
	return nil
}

func anglerads(u, v, w Vec2) float32 {
	var x Vec2
	x.Sub(&u, &w)
	a := sqrt(x.Dot(&x))
	x.Sub(&u, &v)
	b := sqrt(x.Dot(&x))
	x.Sub(&w, &v)
	c := sqrt(x.Dot(&x))

	s := (a + b + c) / 2
	area := sqrt(s * (s - a) * (s - b) * (s - c))

	return (2 * area) / (a * b)
}

func circumcenter(a, b, c Vec2) Vec2 {

	var x Vec2
	x.Add(&a, &b)
	x.DivScalar(2)
	x.Add(&x, &c)
	x.DivScalar(2)
	log.Println(x)
	return x

	//
	// d := 2 * (a[0]*(b[1]-c[1]) + b[0]*(c[1]-a[1]) + c[0]*(a[1]-b[1]))
	// x := (a.Dot(&a)*b[1] - c[1] + b.Dot(&b)*(c[1]-a[1]) + c.Dot(&c)*(a[1]-b[1])) / d
	// y := (a.Dot(&a)*c[0] - b[0] + b.Dot(&b)*(a[0]-c[0]) + c.Dot(&c)*(b[0]-a[0])) / d
	// log.Println(x, y)
	// return Vec2{x, y}

	//
	// a := Mat3{
	// {u[0], u[1], 1},
	// {v[0], v[1], 1},
	// {w[0], w[1], 1},
	// }.Det()
	// b := Mat3{
	// {u[0], u[1], u.Dot(&u)},
	// {v[0], v[1], v.Dot(&v)},
	// {w[0], w[1], w.Dot(&w)},
	// }.Det()
	// log.Println("circumcenter", a, b)
	// return Vec2{a, b}
}

func (st *Store2D) hasboundary(u, v, w Vec2) bool {
	_, a := st.Adjacent(v, u)
	_, b := st.Adjacent(w, v)
	_, c := st.Adjacent(u, w)
	return !a || !b || !c
}

func (st *Store2D) FixAngles() {
	const y = 0.39
	// var err error
LOOP:
	for e0, e1 := range st.m {
		u, v, w := e0[0], e1[0], e1[1]
		if u == GhostVertex || v == GhostVertex || w == GhostVertex {
			continue
		}
		if !st.hasboundary(u, v, w) && anglerads(u, v, w) < y {
			// x := circumcenter(u, v, w)
			// u, v, w, err = st.PointLocation(x)
			// if err == nil {
			// st.InsertVertex(x, u, v, w)
			// goto LOOP
			// }
			st.InsertVertex(circumcenter(u, v, w), u, v, w)
			goto LOOP
			// break
		}
		if !st.hasboundary(v, w, u) && anglerads(v, w, u) < y {
			// x := circumcenter(v, w, u)
			// u, v, w, err = st.PointLocation(x)
			// if err == nil {
			// st.InsertVertex(x, u, v, w)
			// goto LOOP
			// }
			st.InsertVertex(circumcenter(v, w, u), v, w, u)
			goto LOOP
			// break
		}
		if !st.hasboundary(w, u, v) && anglerads(w, u, v) < y {
			// x := circumcenter(w, u, v)
			// u, v, w, err = st.PointLocation(x)
			// if err == nil {
			// st.InsertVertex(x, u, v, w)
			// goto LOOP
			// }
			st.InsertVertex(circumcenter(w, u, v), w, u, v)
			goto LOOP
			// break
		}
	}
}

// func (st *Store2D) InsertVertexAtConflict(u Vec2, v, w, x Vec2) {
// }

// func (st *Store2D) MarkCavity(u Vec2, v, w Vec2, D []Vec2, C []Vec2) {
// x, _ := st.Adjacent(w, v)
// if b := st.v[[3]Vec2{w, v, x}]; b {
// return
// }
// if InCircle(&u, &v, &w, &x) > 0 {
// st.v[[3]Vec2{w, v, x}] = true
// D = append(D, w, v, x)
// st.MarkCavity(u, v, x, D, C)
// st.MarkCavity(u, x, w, D, C)
// } else {
// C = append(C, u, v, w)
// }
// }

// func (st *Store2D) RedistributeList(u Vec2, v, w, x Vec2) {
// for _, y := range st.g[[3]Vec2{v, w, x}] {

// }
// }
