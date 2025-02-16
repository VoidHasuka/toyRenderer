#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>

#include "model.h"
#include "Triangle.h"

Model::Model(const char* filename) :vertNum(0), faceNum(0) {
	std::vector<Vec4f> verts_;
	std::vector<std::vector<Vec3i>> faces_;
	std::vector<Vec3f> norms_;
	std::vector<Vec2f> uv_;

	std::ifstream in;
	in.open(filename, std::ifstream::in);
	if (in.fail())return;
	std::string line;
	while (!in.eof()) {
		std::getline(in, line);
		std::istringstream iss(line.c_str());
		char trash;
		if (!line.compare(0, 2, "v ")) {
			iss >> trash;
			Vec4f v;
			for (int i = 0; i < 3; i++)iss >> v[i];
			v[3] = 1.0f;
			verts_.push_back(v);
		}
		else if (!line.compare(0, 3, "vt ")) {
			iss >> trash >> trash;
			Vec2f uv;
			for (int i = 0; i < 2; i++)iss >> uv[i];
			uv_.push_back(uv);
		}
		else if (!line.compare(0, 3, "vn ")) {
			iss >> trash >> trash;
			Vec3f normal;
			for (int i = 0; i < 3; i++)iss >> normal[i];
			norms_.push_back(normal);
		}
		else if (!line.compare(0, 2, "f ")) {
			std::vector<Vec3i> f;
			Vec3i tmp;
			iss >> trash;
			while (iss >> tmp[0] >> trash >> tmp[1] >> trash >> tmp[2]) {
				for (int i = 0; i < 3; i++)tmp[i]--;
				f.push_back(tmp);
			}
			faces_.push_back(f);
		}
	}
	std::cerr << "# v# " << verts_.size() << " f# " << faces_.size() << " vt# " << uv_.size() << " vn# " << norms_.size() << std::endl;

	TriangleList.resize(faces_.size());
	TriangleList.clear();

	for (int i = 0; i < faces_.size(); i++) {
		Triangle tmp;

		std::vector<Vec3i> face = faces_[i];
		Vec3i face_vertex_1 = face[0];
		Vec3i face_vertex_2 = face[1];
		Vec3i face_vertex_3 = face[2];

		tmp.v[0] = verts_[face_vertex_1.x];
		tmp.v[1] = verts_[face_vertex_2.x];
		tmp.v[2] = verts_[face_vertex_3.x];

		tmp.texCoords[0] = uv_[face_vertex_1.y];
		tmp.texCoords[1] = uv_[face_vertex_2.y];
		tmp.texCoords[2] = uv_[face_vertex_3.y];

		tmp.normal[0] = norms_[face_vertex_1.z];
		tmp.normal[1] = norms_[face_vertex_2.z];
		tmp.normal[2] = norms_[face_vertex_3.z];

		TriangleList.push_back(tmp);
	}

	vertNum = static_cast<int>(verts_.size());
	faceNum = static_cast<int>(faces_.size());
}

Model::~Model() {

}

int Model::nverts() {
	return vertNum;
}

int Model::nfaces() {
	return faceNum;
}
