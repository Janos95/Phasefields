

#include "modica_mortola.hpp"
#include "color_callback.hpp"
#include "phasefield_initialization.hpp"
#include "viewer.hpp"
#include "upload.hpp"
#include "optimization_context.hpp"


using namespace Magnum;
using namespace Corrade;
using namespace Corrade::Containers;

using namespace std::chrono_literals;


class PhasefieldHandle : public ceres::IterationCallback
{
public:
    explicit PhasefieldHandle(const Eigen::VectorXd& U): m_U(U)
    {
    }

    ceres::CallbackReturnType operator()(const ceres::IterationSummary&) override
    {
        std::lock_guard l(m_mutex);
        m_UCopy = m_U;
        current = true;
        return ceres::CallbackReturnType::SOLVER_CONTINUE;
    }

    bool updateIfCurrent(Eigen::VectorXd& U) {
        std::lock_guard l(m_mutex);
        if(current){
            U = m_UCopy;
            current = false;
            return true;
        }
        return false;
    }

private:

    bool current = false;
    Corrade::Containers::Reference<const Eigen::VectorXd> m_U;
    Eigen::VectorXd m_UCopy;

    mutable std::mutex m_mutex;
};

int main(int argc, char** argv) {
    
    //ScopedTimer ti("total");

    //Utility::Configuration conf(CONF_PATH, Utility::Configuration::Flag::ReadOnly);

    //auto epsilon = conf.value<double>("epsilon");
    //auto a = conf.value<double>("a");
    //auto b = conf.value<double>("b");
    //auto wMM = conf.value<double>("modica_mortola_weight");
    //auto wArea = conf.value<double>("area_weight");
    //auto wConnectedness = conf.value<double>("connectedness_weight");
    //auto enforceConnectedness = conf.value<bool>("enforce_connectedness");
    //auto postConnect = conf.value<bool>("post_connect");
    //auto iterations = conf.value<int>("iterations");
    //auto postConnectIterations = conf.value<int>("post_connect_iterations");
    //auto inputPath = conf.value<std::string>("input_path");
    //auto outputPath1 = conf.value<std::string>("output_path1");
    //auto outputPath2 = conf.value<std::string>("output_path2");

    //auto interfaceLoss = std::make_unique<ceres::TrivialLoss>();
    //auto potentialLoss = std::make_unique<ceres::TrivialLoss>();
    //auto areaLoss = std::make_unique<ceres::TrivialLoss>();
    //auto connectednessLoss1 = std::make_unique<ceres::TrivialLoss>();
    //auto connectednessLoss2 = std::make_unique<ceres::TrivialLoss>();

    //auto meshdata = loadMeshData(inputPath);

    //auto [V, F] = toEigen(meshdata);

    //fmt::print("Solving Problem with {} vertices\n", V.rows());

    //Eigen::VectorXd U(V.rows());
    //initRandomNormal(U);

    //auto* cost = new SumProblem(U.size());
    //cost->push_back(std::make_unique<InterfaceEnergy>(V,F,epsilon), std::move(interfaceLoss), wMM);
    //cost->push_back(std::make_unique<PotentialEnergy>(V,F,epsilon), std::move(potentialLoss), wMM);
    //cost->push_back(std::make_unique<AreaRegularizer>(V,F), std::move(areaLoss), wArea);
    //if(enforceConnectedness){
    //    cost->push_back(std::make_unique<ConnectednessConstraint>(V, F, epsilon, a, b), std::move(connectednessLoss1), wConnectedness);
    //    cost->push_back(std::make_unique<ConnectednessConstraint>(V, F, epsilon, -b, -a), std::move(connectednessLoss2), wConnectedness);
    //}

    //PhasefieldHandle phasefieldHandle(U);

    //// Run the solver!
    //ceres::GradientProblemSolver::Options options;
    //options.minimizer_progress_to_stdout = true;
    //options.max_num_iterations = iterations;
    //options.update_state_every_iteration = true;
    //options.callbacks.push_back(&phasefieldHandle);
    //ceres::GradientProblemSolver::Summary summary;
    //std::thread t([&] {
    //                  ceres::GradientProblem problem(cost);
    //                  ceres::Solve(options, problem, U.data(), &summary);
    //              });

    //Viewer viewer(argc, argv);
    //Scene scene;
    //viewer.setScene(scene);

    //auto* object = scene.addObject("mesh", upload(std::move(meshdata), CompileFlag::GenerateSmoothNormals));
    //auto& generated = object->meshData;

    //scene.addNode("mesh", "mesh", ShaderType::Phong);
    //viewer.callbacks.emplace_back([&]
    //    (auto& v) mutable
    //    {
    //        Eigen::VectorXd C;
    //        if(phasefieldHandle.updateIfCurrent(C)){
    //            auto* obj = scene.getObject("mesh");
    //            CORRADE_INTERNAL_ASSERT(obj);
    //            GL::Buffer& vertices = obj->vertices;
    //            auto data = vertices.map(
    //                    0,
    //                    vertices.size(),
    //                    GL::Buffer::MapFlag::Write);

    //            CORRADE_CONSTEXPR_ASSERT(data, "could not map vertex data");
    //            auto colorView = generated.mutableAttribute<Color4>(Trade::MeshAttribute::Color);
    //            std::transform(C.begin(), C.end(), colorView.begin(), JetColorMap{});

    //            Utility::copy(generated.vertexData(), data);
    //            vertices.unmap();
    //            v.redraw();
    //        }
    //    }
    //);

    Viewer viewer(argc, argv);
    OptimizationContext optimizationContext;
    viewer.tickCallbacks.emplace_back([&](auto& s){ optimizationContext.updateScene(s); });
    viewer.menuCallbacks.emplace_back([&](auto& c){ optimizationContext.showMenu(c); });
    viewer.exec();
}
