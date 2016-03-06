package tds

import "testing"

func TestOrient2D(t *testing.T) {
	v0, v1, v2 := Vec2{0, 0}, Vec2{1, 0}, Vec2{1, 1}
	if x := Orient2D(v0, v1, v2); x != 1 {
		t.Errorf("want 1, have %v", x)
	}
	if x := Orient2D(v0, v2, v1); x != -1 {
		t.Errorf("want -1, have %v", x)
	}
	if x := Orient2D(v0, v1, v0); x != 0 {
		t.Errorf("want 0, have %v", x)
	}
}

func TestOrient3D(t *testing.T) {
	v0, v1, v2, v3 := Vec3{0, 1, 0}, Vec3{0, 0, 1}, Vec3{1, 0, 0}, Vec3{0, 0, -1}
	if x := Orient3D(v0, v1, v2, v3); x != 2 {
		t.Errorf("want 2, have %v", x)
	}
	if x := Orient3D(v0, v2, v1, v3); x != -2 {
		t.Errorf("want -2, have %v", x)
	}
	if x := Orient3D(v0, v1, v0, v3); x != 0 {
		t.Errorf("want 0, have %v", x)
	}
}

func TestInCircle(t *testing.T) {
	v0, v1, v2 := &Vec2{0, 0}, &Vec2{1, 0}, &Vec2{1, 1}
	if x := Orient2D(*v0, *v1, *v2); x <= 0 {
		t.Fatal("vectors not oriented counter-clockwise")
	}
	if x := InCircle(v0, v1, v2, &Vec2{0.5, 0.5}); x != 0.5 {
		t.Errorf("want 0.5, have %v", x)
	}
	if x := InCircle(v0, v1, v2, &Vec2{-1, 0}); x != -2 {
		t.Errorf("want -2, have %v", x)
	}
	if x := InCircle(v0, v1, v2, &Vec2{0, 1}); x != 0 {
		t.Errorf("want 0, have %v", x)
	}
}

func TestInSphere(t *testing.T) {
	v0, v1, v2, v3 := Vec3{0, 1, 0}, Vec3{0, 0, 1}, Vec3{1, 0, 0}, Vec3{0, 0, -1}
	if x := Orient3D(v0, v1, v2, v3); x <= 0 {
		t.Fatal("vectors not oriented counter-clockwise")
	}
	if x := InSphere(v0, v1, v2, v3, Vec3{1, 1, 1}); x != -4 {
		t.Errorf("want -4, have %v", x)
	}
	if x := InSphere(v0, v1, v2, v3, Vec3{0, 0, 0}); x != 2 {
		t.Errorf("want 2, have %v", x)
	}
	if x := InSphere(v0, v1, v2, v3, Vec3{0, -1, 0}); x != 0 {
		t.Errorf("want 0, have %v", x)
	}
}

func TestStore2D(t *testing.T) {
	st := NewStore2D()
	v0, v1, v2 := Vec2{0, 0}, Vec2{1, 0}, Vec2{1, 1}
	if err := st.AddTriangle(v0, v1, v2); err != nil {
		t.Error(err)
	}
	if err := st.DeleteTriangle(v0, v1, v2); err != nil {
		t.Error(err)
	}
	if err := st.AddTriangle(v0, v1, v2); err != nil {
		t.Error(err)
	}
}
