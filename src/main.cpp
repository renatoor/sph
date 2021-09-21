#include <Corrade/Containers/Optional.h>
#include <Magnum/GL/DefaultFramebuffer.h>
#include <Magnum/GL/Renderer.h>
#include <Magnum/Platform/GlfwApplication.h>
#include <Magnum/Timeline.h>
#include <Magnum/ImGuiIntegration/Context.hpp>

#include "Objects/ArcBall.h"
#include "Objects/ArcBallCamera.h"
#include "Objects/WireframeObjects.h"
#include "Objects/ParticleGroup.h"
#include "SPH.h"

namespace Magnum { namespace Examples {

using Object3D = SceneGraph::Object<SceneGraph::MatrixTransformation3D>;
using Scene3D = SceneGraph::Scene<SceneGraph::MatrixTransformation3D>;

using namespace Math::Literals;

namespace {
    constexpr Float ParticleRadius = 0.0457f;
}

int count = 0;

class MainApplication: public Platform::Application {
    public:
        explicit MainApplication(const Arguments& arguments);

    private:
        void drawEvent() override;
        void viewportEvent(ViewportEvent& event) override;
        void keyPressEvent(KeyEvent& event) override;
        void mousePressEvent(MouseEvent& event) override;
        void mouseReleaseEvent(MouseEvent& event) override;
        void mouseMoveEvent(MouseMoveEvent& event) override;
        void mouseScrollEvent(MouseScrollEvent& event) override;

        Containers::Pointer<Scene3D> _scene;
        Containers::Pointer<SceneGraph::DrawableGroup3D> _drawables;

        Containers::Optional<ArcBallCamera> _arcballCamera;

        Containers::Pointer<WireframeBox> _box;

        Containers::Pointer<ParticleGroup> _drawableParticles;
        
        std::vector<Vector3> _particlePositions;
        //std::vector<Vector3> _particlePositions = {
        //    Vector3(0.0f, 0.7f, 0.0f),
        //    Vector3(0.0f, 0.9f, 0.0f),
        //};

        SPH *_sph;

        Timeline _timeline;

        ImGuiIntegration::Context _imgui{NoCreate};

        bool _temperature = true;
        bool _density = true;
        bool _velocity = true;
        bool _acceleration = true;
        bool _position = true;
};

MainApplication::MainApplication(const Arguments& arguments) :
    Platform::Application{arguments, NoCreate}
{
    /* Setup window */
    {
        const Vector2 dpiScaling = this->dpiScaling({});

        Configuration conf;

        conf.setTitle("SPH")
            .setSize(conf.size(), dpiScaling)
            .setWindowFlags(Configuration::WindowFlag::Resizable);

        GLConfiguration glConf;

        glConf.setSampleCount(dpiScaling.max() < 2.0f ? 8 : 2);

        if (!tryCreate(conf, glConf)) {
            create(conf, glConf.setSampleCount(0));
        }
    }

    _sph = new SPH(ParticleRadius, 5000);
    /*
    float min = 0.25f;
    float max = 0.75f;
    float med = 0.0f;
    float inc = 0.04f;

    for (float x = 0.75f - 0.04f; x >= 1.0f - max; x -= inc) {
        for (float y = 0.1f; y <= 0.8f; y += inc) {
            for (float z = 0.4; z < 0.6f; z += inc) {
                float jitter = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
                jitter /= 1000.0f;
                _sph->addParticle(Vector3(
                            Math::clamp(x + jitter, 0.0f, 1.0f),
                            y,
                            Math::clamp(z + jitter, 0.0f, 1.0f)
                            ), Vector3(0.0f));
            }
        }
    }
    */

    _scene.reset(new Scene3D{});
    _drawables.reset(new SceneGraph::DrawableGroup3D{});
 
    /* Set up the camera */
    {
        /* Setup the arcball after the camera objects */
        const Vector3 eye = Vector3(0.5f, 0.5f, -2.0f);
        const Vector3 center{0.5f, 0.5f, 0.5f};
        const Vector3 up = Vector3::yAxis();

        _arcballCamera.emplace(_scene.get(), eye, center, up, 45.0_degf,
            windowSize(), framebufferSize());
    }

    _box.reset(new WireframeBox(_scene.get(), _drawables.get()));
    _box->transform(Matrix4::scaling(Vector3{ 0.5, 0.5, 0.5 }) * Matrix4::translation(Vector3(1)));
    _box->setColor(Color3(1, 1, 0));

    //_drawableParticles.reset(new ParticleGroup{_sph->particlePositions(), ParticleRadius});
    _drawableParticles.reset(new ParticleGroup{_sph->particleVertex(), ParticleRadius});
    _drawableParticles->setDirty();

    _drawableParticles->setColorMode(ParticleSphereShader::ColorMode(0));

    GL::Renderer::enable(GL::Renderer::Feature::DepthTest);
    GL::Renderer::enable(GL::Renderer::Feature::FaceCulling);
    GL::Renderer::enable(GL::Renderer::Feature::ProgramPointSize);

    _imgui = ImGuiIntegration::Context(Vector2{windowSize()}/dpiScaling(),
            windowSize(), framebufferSize());

    /* Set up proper blending to be used by ImGui. There's a great chance
       you'll need this exact behavior for the rest of your scene. If not, set
       this only for the drawFrame() call. */
    GL::Renderer::setBlendEquation(GL::Renderer::BlendEquation::Add,
            GL::Renderer::BlendEquation::Add);
    GL::Renderer::setBlendFunction(GL::Renderer::BlendFunction::SourceAlpha,
            GL::Renderer::BlendFunction::OneMinusSourceAlpha);

    setSwapInterval(0);

    _timeline.start();
}

void MainApplication::drawEvent() {
    GL::defaultFramebuffer.clear(GL::FramebufferClear::Color);

    Float deltaTs = _timeline.previousFrameDuration();

    if (count < 1500) {
        float jitter = static_cast<float>(rand()) / static_cast<float>(RAND_MAX) / 5.f;
        float jitter2 = static_cast<float>(rand()) / static_cast<float>(RAND_MAX) / 5.f;
        //_sph->addParticle(Vector3(1.0f, 0.8f + jitter, 0.5f + jitter), Vector3(3.0f + jitter2, 1.0f, jitter));
        _sph->addParticle(Vector3(0.96f - jitter, 0.04f, 0.96f - jitter2), Vector3(-0.8f, 0.0f, -0.8f), 25.0f);
        _sph->addParticle(Vector3(0.04f + jitter, 0.04f, 0.04f + jitter2), Vector3(1.0f, 0.0f, 1.0f), 35.0f);
        count++;
    }

    GL::defaultFramebuffer.clear(GL::FramebufferClear::Color | GL::FramebufferClear::Depth);

    _sph->step(deltaTs);

    _drawableParticles->setDirty();
    /* Draw particles */
    _drawableParticles->draw(_arcballCamera->camera(), framebufferSize());

    /* Call arcball update in every frame. This will do nothing if the camera
       has not been changed. Otherwise, camera transformation will be
       propagated into the camera objects. */

    //bool camChanged = _arcballCamera->update();
    _arcballCamera->update();
    _arcballCamera->draw(*_drawables);

    _imgui.newFrame();

    ImGui::Checkbox("Temperatura", &_temperature);

    if (_temperature) {
        _sph->showTemperature();
        _density = false;
        _velocity = false;
        _acceleration = false;
        _position = false;
    }

    ImGui::Checkbox("Densidade", &_density);

    if (_density) {
        _sph->showDensity();
        _temperature = false;
        _velocity = false;
        _acceleration = false;
        _position = false;
    }

    ImGui::Checkbox("Velocidade", &_velocity);

    if (_velocity) {
        _sph->showVelocity();
        _temperature = false;
        _density = false;
        _acceleration = false;
        _position = false;
    }

    ImGui::Checkbox("Aceleração", &_acceleration);

    if (_acceleration) {
        _sph->showAcceleration();
        _temperature = false;
        _density = false;
        _velocity = false;
        _position = false;
    }

    ImGui::Checkbox("Posição", &_position);

    if (_position) {
        _sph->showPosition();
        _temperature = false;
        _density = false;
        _velocity = false;
        _acceleration = false;
    }

    _imgui.updateApplicationCursor(*this);

    GL::Renderer::enable(GL::Renderer::Feature::Blending);
    GL::Renderer::enable(GL::Renderer::Feature::ScissorTest);
    GL::Renderer::disable(GL::Renderer::Feature::FaceCulling);
    GL::Renderer::disable(GL::Renderer::Feature::DepthTest);

    _imgui.drawFrame();

    GL::Renderer::enable(GL::Renderer::Feature::DepthTest);
    GL::Renderer::enable(GL::Renderer::Feature::FaceCulling);
    GL::Renderer::disable(GL::Renderer::Feature::ScissorTest);
    GL::Renderer::disable(GL::Renderer::Feature::Blending);


    swapBuffers();

    //if(camChanged) redraw();
    redraw();
    _timeline.nextFrame();
}

void MainApplication::viewportEvent(ViewportEvent& event) {
    GL::defaultFramebuffer.setViewport({{}, event.framebufferSize()});

    _arcballCamera->reshape(event.windowSize(), event.framebufferSize());
}

void MainApplication::keyPressEvent(KeyEvent& event) {
    switch(event.key()) {
        case KeyEvent::Key::L:
            if(_arcballCamera->lagging() > 0.0f) {
                Debug{} << "Lagging disabled";
                _arcballCamera->setLagging(0.0f);
            } else {
                Debug{} << "Lagging enabled";
                _arcballCamera->setLagging(0.85f);
            }
            break;
        case KeyEvent::Key::R:
            _arcballCamera->reset();
            break;

        default: return;
    }

    event.setAccepted();
    redraw(); /* camera has changed, redraw! */
}

void MainApplication::mousePressEvent(MouseEvent& event) {
    _imgui.handleMousePressEvent(event);
    _arcballCamera->initTransformation(event.position());

    event.setAccepted();
    redraw(); /* camera has changed, redraw! */
}

void MainApplication::mouseReleaseEvent(MouseEvent& event) {
    _imgui.handleMouseReleaseEvent(event);
}

void MainApplication::mouseMoveEvent(MouseMoveEvent& event) {
    _imgui.handleMouseMoveEvent(event);
    if(!event.buttons()) return;

    if(event.modifiers() & MouseMoveEvent::Modifier::Shift)
        _arcballCamera->translate(event.position());
    else _arcballCamera->rotate(event.position());

    event.setAccepted();
    redraw(); /* camera has changed, redraw! */
}

void MainApplication::mouseScrollEvent(MouseScrollEvent& event) {
    const Float delta = event.offset().y();
    if(Math::abs(delta) < 1.0e-2f) return;

    _arcballCamera->zoom(delta);

    event.setAccepted();
    redraw(); /* camera has changed, redraw! */
}

}}

MAGNUM_APPLICATION_MAIN(Magnum::Examples::MainApplication)
