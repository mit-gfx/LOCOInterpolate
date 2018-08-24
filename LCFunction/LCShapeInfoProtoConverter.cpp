#include <Eigen/Dense>
#include <fstream>
#include <sstream>
#include <iostream>
#include "LCShapeInfoProtoConverter.h"
#include "LCFunctionValue.h"
#include "LCError.h"

//-------------------------------------------------------------------------------

/*
LCError OCTetCADShapeInfoConverter::outputToProto(std::string filename, OCTetCADShapeInfo* shape)
{
	LCError err;
	//std::cout << "before creating proto" << std::endl;
	//system("pause");
	proto::TetMesh  shapeProto;
	LCErrorReturn(convertToProto(shape, &shapeProto));
	LCErrorReturn(saveToFile(shapeProto, filename.c_str(), false));
	//std::cout << "after creating proto" << std::endl;
	//system("pause");

	return err;
}
*/

/*
LCError OCTetCADShapeInfoConverter::loadFromProto(std::string filename, OCTetCADShapeInfo* shape)
{
	//std::cout << "before creating proto" << std::endl;
	//system("pause");
	LCError err;
	proto::TetMesh shapeProto;
	LCErrorReturn(loadFromFile(filename.c_str(), false, &shapeProto));
	LCErrorReturn(convertToCplusplus(shapeProto, shape));
	//std::cout << "after creating proto" << std::endl;
	//system("pause");
	return err;

}
*/

/*
void convertConditionValueToCpp(const proto::BoundaryConditionValue &proto, OCBoundaryValue * cpp)
{
	switch (proto.boundarytype())
	{
	case proto::BoundaryType::FIXED_BOUNDARY:
		cpp->type = OCBoundaryValue::FIXED_BOUNDARY;
		break;
	case proto::BoundaryType::BOUNDARY_FORCE:
		cpp->type = OCBoundaryValue::BOUNDARY_FORCE;
		break;
	case proto::BoundaryType::BOUNDARY_HEAT:
		cpp->type = OCBoundaryValue::BOUNDARY_HEAT;
		break;
	default:
		break;
	}
	proto::Vector3d protoVer = proto.force();
	cpp->force =  Eigen::Vector3d(protoVer.x(), protoVer.y(), protoVer.z());
	cpp->heat = proto.heatvalue();
}

void convertConditionValueToProto(const OCBoundaryValue &cpp, proto::BoundaryConditionValue * proto)
{
	switch (cpp.type)
	{
	case OCBoundaryValue::FIXED_BOUNDARY:
		proto->set_boundarytype(proto::BoundaryType::FIXED_BOUNDARY);
		break;
	case OCBoundaryValue::BOUNDARY_FORCE:
		proto->set_boundarytype(proto::BoundaryType::BOUNDARY_FORCE);
		break;
	case OCBoundaryValue::BOUNDARY_HEAT:
		proto->set_boundarytype(proto::BoundaryType::BOUNDARY_HEAT);
		break;
	default:
		break;
	}

	proto::Vector3d * force = new proto::Vector3d();
	force->set_x(cpp.force(0));
	force->set_y(cpp.force(1));
	force->set_z(cpp.force(2));
	proto->set_allocated_force(force);
	proto->set_heatvalue(cpp.heat);
}


LCError OCTetCADShapeInfoConverter::convertToProto(OCTetCADShapeInfo* shapeInfo, proto::TetMesh* result)
{
	LCError err;
	//std::cout << "converting tets" << std::endl;
	for (auto tet : shapeInfo->tetMesh.tets)
	{
		proto::Vector4i * vec = result->add_tets();
		vec->set_x(tet(0));
		vec->set_y(tet(1));
		vec->set_z(tet(2));
		vec->set_w(tet(3));
	}
	//std::cout << "converting vertex" << std::endl;
	for (auto vertex : shapeInfo->tetMesh.vertices)
	{
		proto::Vector3d * vec = result->add_vertices();
		vec->set_x(vertex(0));
		vec->set_y(vertex(1));
		vec->set_z(vertex(2));
	}
	//std::cout << "converting boundary" << std::endl;

	for (auto boudary : shapeInfo->boundaryConditions)
	{
		proto::BoundaryCondition * protoBoundary = result->add_boundaryconditions();
		convertConditionValueToProto(boudary->value, protoBoundary->mutable_value());
		for (auto vertex : boudary->vertices)
		{
			protoBoundary->add_vertices(vertex);
		}

	}
	//std::cout << "converting physics" << std::endl;

	for (auto physics : shapeInfo->precomputedPhysics)
	{
		proto::PrecomputedPhysics * protoPhysics = result->add_precomputedphysics();
		protoPhysics->set_name(physics->name);
		for (auto vertex : physics->values)
		{
			protoPhysics->add_values(vertex);
		}
	}
	
	// save the control points
	//std::cout << "converting contol points" << std::endl;

	proto::ControlPoints * protoControlPoints = result->mutable_controlpoints();
	for (auto point : shapeInfo->controlPoints.points)
	{
		proto::ReferencePoint * refPoints = protoControlPoints->add_point();
		refPoints->set_id(point.first);
		proto::Vector3d * vec = refPoints->mutable_pos();
		vec->set_x(point.second(0));
		vec->set_y(point.second(1));
		vec->set_z(point.second(2));
	}
	//std::cout << "converting links" << std::endl;

	for (auto link : shapeInfo->controlPoints.links)
	{
		proto::Link * protoLink = protoControlPoints->add_links();
		protoLink->set_first(link.first);
		protoLink->set_second(link.second);
		for (auto midPoint : link.midPoints)
		{
			proto::Vector3d *vec =  protoLink->add_midpoints();
			vec->set_x(midPoint(0));
			vec->set_y(midPoint(1));
			vec->set_z(midPoint(2));
		}
	}
	//std::cout << "converting boudnary" << std::endl;

	for (auto boundaryControl : shapeInfo->controlPoints.boundaryControls)
	{
		proto::BoundaryControl * protoBoundaryControl = protoControlPoints->add_boundarycontrols();
		for (auto loop : boundaryControl->controlLoops)
		{
			proto::ControlLoop * protoLoop = protoBoundaryControl->add_controlloops();
			for (auto loopId : loop)
			{
				protoLoop->add_loop(loopId);
			}			
		}
		convertConditionValueToProto(boundaryControl->value, protoBoundaryControl->mutable_value());
	}
	//std::cout << "done" << std::endl;

	
	return err;
}


LCError OCTetCADShapeInfoConverter::convertToCplusplus(const proto::TetMesh & tetMesh,
	OCTetCADShapeInfo* result)
{
	//std::cout << "--------read shape = " << result->getFilename() << std::endl;

	for (int i = 0; i < tetMesh.vertices_size(); i++)
	{
		proto::Vector3d protoVer = tetMesh.vertices(i);
		Eigen::Vector3d vertex(protoVer.x(), protoVer.y(), protoVer.z());
		result->tetMesh.vertices.push_back(vertex);
	}
	for (int i = 0; i < tetMesh.tets_size(); i++)
	{
		proto::Vector4i protoVer = tetMesh.tets(i);
		Eigen::Vector4i vertex(protoVer.x(), protoVer.y(), protoVer.z(), protoVer.w());
		result->tetMesh.tets.push_back(vertex);
	}

	for (int i = 0; i < tetMesh.boundaryconditions_size(); i++)
	{
		OCBoundaryCondition * boundaryCondition = new OCBoundaryCondition();
		auto protoBoundary = tetMesh.boundaryconditions(i);
		if (protoBoundary.has_value())
		{
			convertConditionValueToCpp(protoBoundary.value(), &boundaryCondition->value);
		}
		else // old data structure
		{
			proto::Vector3d force = protoBoundary.force();
			boundaryCondition->value.force = Eigen::Vector3d(force.x(), force.y(), force.z());
			if (protoBoundary.isfixed())
			{
				boundaryCondition->value.type = OCBoundaryValue::OCBoundaryType::FIXED_BOUNDARY;
			}
			else
			{
				boundaryCondition->value.type = OCBoundaryValue::OCBoundaryType::BOUNDARY_FORCE;
			}
		}
		
		for (int j = 0; j < protoBoundary.vertices_size(); j++)
		{
			boundaryCondition->vertices.push_back(protoBoundary.vertices(j));
		}
		result->boundaryConditions.push_back(boundaryCondition);
	}
	for (int i = 0; i < tetMesh.precomputedphysics_size(); i++)
	{
		OCPrecomputedPhysics *precomputedPhysics = new OCPrecomputedPhysics();
		auto protoPhysics = tetMesh.precomputedphysics(i);

		precomputedPhysics->name = protoPhysics.name();
		for (int j = 0; j < protoPhysics.values_size(); j++)
		{
			precomputedPhysics->values.push_back(protoPhysics.values(j));
		}

		result->precomputedPhysics.push_back(precomputedPhysics);
	}


	// Load control Points
	if (!tetMesh.has_controlpoints())
	{
		std::cout << "mesh does not have control points" << std::endl;
		return LCError();
	}

	//std::cout << "has the control points" << std::endl;
	for (int i = 0; i < tetMesh.controlpoints().point_size(); i++)
	{
		std::string pointId = tetMesh.controlpoints().point(i).id();
		proto::Vector3d protoVer = tetMesh.controlpoints().point(i).pos();
		Eigen::Vector3d pointVec(protoVer.x(), protoVer.y(), protoVer.z());
		result->controlPoints.points[pointId] = pointVec;
	}
	for (int i = 0; i < tetMesh.controlpoints().links_size(); i++)
	{
		OCLink link;
		link.first = tetMesh.controlpoints().links(i).first();
		link.second = tetMesh.controlpoints().links(i).second();
		for (int j = 0; j < tetMesh.controlpoints().links(i).midpoints_size(); j++)
		{
			proto::Vector3d protoVer = tetMesh.controlpoints().links(i).midpoints(j);
			Eigen::Vector3d vertex(protoVer.x(), protoVer.y(), protoVer.z());
			link.midPoints.push_back(vertex);
		}
		result->controlPoints.links.push_back(link);
	}
	for (int i = 0; i < tetMesh.controlpoints().boundarycontrols_size(); i++)
	{
		OCBoundaryControl * boundaryControl  = new OCBoundaryControl;
		auto loops  = tetMesh.controlpoints().boundarycontrols(i).controlloops();
		for (int j = 0; j < tetMesh.controlpoints().boundarycontrols(i).controlloops_size(); j++)
		{
			std::vector<std::string > loop;
			auto protoLoop = tetMesh.controlpoints().boundarycontrols(i).controlloops(j);
			for (int k = 0; k < protoLoop.loop_size(); k++)
			{
				loop.push_back(protoLoop.loop(k));
			}
			boundaryControl->controlLoops.push_back(loop);
		}		
		convertConditionValueToCpp(tetMesh.controlpoints().boundarycontrols(i).value(), &boundaryControl->value);

		result->controlPoints.boundaryControls.push_back(boundaryControl);		

	}
	//result->controlPoints.log();
	//std::cout << "-------------------------------------------" << std::endl;

	return LCError();
}
*/
