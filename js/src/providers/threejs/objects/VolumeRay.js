// jshint maxstatements:false,maxcomplexity:false

const THREE = require('three');
const _ = require('../../../lodash');
const colorMapHelper = require('../../../core/lib/helpers/colorMap');
const { closestPowOfTwo } = require('../helpers/Fn');
const { typedArrayToThree } = require('../helpers/Fn');
const { areAllChangesResolve } = require('../helpers/Fn');
const { commonUpdate } = require('../helpers/Fn');

/**
 * Loader strategy to handle Volume object
 * @method Volume
 * @memberof K3D.Providers.ThreeJS.Objects
 * @param {Object} config all configurations params from JSON
 * @param {K3D}
 * @return {Object} 3D object ready to render
 */
module.exports = {
    create(config, K3D) {
        config.samples = config.samples || 512.0;
        config.alpha_coef = typeof (config.alpha_coef) !== 'undefined' ? config.alpha_coef : 50.0;
        config.gradient_step = config.gradient_step || 0.005;
        config.shadow = config.shadow || 'off';
        config.interpolation = typeof (config.interpolation) !== 'undefined' ? config.interpolation : true;
        config.shadow_delay = config.shadow_delay || 500;
        config.shadow_res = closestPowOfTwo(config.shadow_res || 128);

        config.ray_samples_count = config.ray_samples_count || 16;
        config.focal_plane = config.focal_plane || 512.0;
        config.focal_length = typeof (config.focal_length) !== 'undefined' ? config.focal_length : 0.0;

        const randomMul = typeof (window.randomMul) !== 'undefined' ? window.randomMul : 255.0;
        const gl = K3D.getWorld().renderer.getContext();
        const geometry = new THREE.BoxBufferGeometry(1, 1, 1);
        const modelMatrix = new THREE.Matrix4();
        const translation = new THREE.Vector3();
        const rotation = new THREE.Quaternion();
        const scale = new THREE.Vector3();
        const maxTextureSize = gl.getParameter(gl.MAX_TEXTURE_SIZE);
        let lightMapSize = config.shadow_res;
        const colorMap = (config.color_map && config.color_map.data) || null;
        let opacityFunction = (config.opacity_function && config.opacity_function.data) || null;
        const colorRange = config.color_range;
        const { samples } = config;
        let sceneRTT;
        let cameraRTT;
        let quadRTT;
        let textureRTT;
        let jitterTexture;
        let listenersId;
        let timeoutId;
        let lastShadowMapUpdated = 0;

        lightMapSize = (lightMapSize > 512 ? 512 : lightMapSize);
        const lightMapRenderTargetSize = closestPowOfTwo(Math.sqrt(lightMapSize * lightMapSize * lightMapSize));

        if (lightMapRenderTargetSize > maxTextureSize) {
            throw new Error(`To big light map size. gl.MAX_TEXTURE_SIZE=${maxTextureSize}`);
        }

        if (opacityFunction === null) {
            opacityFunction = [colorMap[0], 0.0, colorMap[colorMap.length - 4], 1.0];
        }

        modelMatrix.set.apply(modelMatrix, config.model_matrix.data);
        modelMatrix.decompose(translation, rotation, scale);

        console.log('Creating data texture');
        const dataTexture = new THREE.DataTexture3D(
            config.volume.data,
            config.volume.shape[2],
            config.volume.shape[1],
            config.volume.shape[0],
        );
        const childTextureLowerBits = new THREE.DataTexture3D(
            config.children_lower_bits.data,
            config.children_lower_bits.shape[2],
            config.children_lower_bits.shape[1],
            config.children_lower_bits.shape[0],
        );
        const childTextureUpperBits = new THREE.DataTexture3D(
            config.children_upper_bits.data,
            config.children_upper_bits.shape[2],
            config.children_upper_bits.shape[1],
            config.children_upper_bits.shape[0],
        );

        const textures = [dataTexture, childTextureLowerBits, childTextureUpperBits];
        textures.forEach(
            (texture) => {
                texture.format = THREE.RedFormat;
                texture.type = typedArrayToThree(config.volume.data.constructor);

                texture.generateMipmaps = false;

                if (config.interpolation) {
                    texture.minFilter = THREE.LinearFilter;
                    texture.magFilter = THREE.LinearFilter;
                } else {
                    texture.minFilter = THREE.NearestFilter;
                    texture.magFilter = THREE.NearestFilter;
                }

                texture.wrapS = THREE.ClampToEdgeWrapping;
                texture.wrapT = THREE.ClampToEdgeWrapping;
                texture.wrapR = THREE.ClampToEdgeWrapping;
                texture.needsUpdate = true;
            },
        );

        jitterTexture = new THREE.DataTexture(
            new Uint8Array(_.range(64 * 64).map(() => randomMul * Math.random())),
            64, 64, THREE.RedFormat, THREE.UnsignedByteType,
        );
        jitterTexture.minFilter = THREE.LinearFilter;
        jitterTexture.magFilter = THREE.LinearFilter;
        jitterTexture.wrapS = THREE.MirroredRepeatWrapping;
        jitterTexture.wrapT = THREE.MirroredRepeatWrapping;
        jitterTexture.generateMipmaps = false;
        jitterTexture.needsUpdate = true;

        const canvas = colorMapHelper.createCanvasGradient(colorMap, 1024, opacityFunction);
        const colormap = new THREE.CanvasTexture(canvas, THREE.UVMapping, THREE.ClampToEdgeWrapping,
            THREE.ClampToEdgeWrapping, THREE.NearestFilter, THREE.NearestFilter);
        colormap.needsUpdate = true;

        if (config.shadow !== 'off') {
            textureRTT = new THREE.WebGLRenderTarget(lightMapRenderTargetSize, lightMapRenderTargetSize, {
                minFilter: THREE.LinearFilter,
                magFilter: THREE.LinearFilter,
                format: THREE.RedFormat,
                type: THREE.UnsignedByteType,
                generateMipmaps: false,
                stencilBuffer: false,
                depthBuffer: false,
            });
        }

        const uniforms = {
            lightMapSize: { value: new THREE.Vector3(lightMapSize, lightMapSize, lightMapSize) },
            volumeMapSize: {
                value: new THREE.Vector3(config.volume.shape[2],
                    config.volume.shape[1],
                    config.volume.shape[0]),
            },
            lightMapRenderTargetSize: { value: new THREE.Vector2(lightMapRenderTargetSize, lightMapRenderTargetSize) },
            low: { value: colorRange[0] },
            high: { value: colorRange[1] },
            samples: { value: samples },
            alpha_coef: { value: config.alpha_coef },
            gradient_step: { value: config.gradient_step },
            translation: { value: translation },
            rotation: { value: rotation },
            shadowTexture: { type: 't', value: (textureRTT ? textureRTT.texture : null) },
            focal_length: { value: config.focal_length },
            focal_plane: { value: config.focal_plane },
            scale: { value: scale },
            volumeTexture: { type: 't', value: dataTexture },
            childTextureLowerBits: { type: 't', value: childTextureLowerBits },
            childTextureUpperBits: { type: 't', value: childTextureUpperBits },
            Npack: { type: 't', value: config.Npack },
            colormap: { type: 't', value: colormap },
            jitterTexture: { type: 't', value: jitterTexture },
        };

        const material = new THREE.ShaderMaterial({
            uniforms: _.merge(
                uniforms,
                THREE.UniformsLib.lights,
            ),
            defines: {
                USE_SHADOW: (config.shadow !== 'off' ? 1 : 0),
                RAY_SAMPLES_COUNT: config.focal_length !== 0.0 ? config.ray_samples_count : 0,
            },
            vertexShader: require('./shaders/Volume.vertex.glsl'),
            fragmentShader: require('./shaders/VolumeRay.fragment.glsl'),
            side: THREE.BackSide,
            depthTest: false,
            depthWrite: false,
            lights: true,
            clipping: true,
            transparent: true,
        });

        geometry.computeBoundingSphere();
        geometry.computeBoundingBox();

        const object = new THREE.Mesh(geometry, material);
        object.applyMatrix4(modelMatrix);
        object.updateMatrixWorld();

        /*
            Light Map support
         */
        if (config.shadow !== 'off') {
            sceneRTT = new THREE.Scene();
            quadRTT = new THREE.Mesh(
                new THREE.PlaneBufferGeometry(lightMapRenderTargetSize, lightMapRenderTargetSize),
                new THREE.ShaderMaterial({
                    uniforms: _.merge(
                        uniforms,
                        THREE.UniformsLib.lights,
                        {
                            lightDirection: { type: 'v3', value: new THREE.Vector3() },
                        },
                    ),
                    defines: {
                        USE_MAP: 1,
                    },
                    vertexShader: require('./shaders/Volume.lightmap.vertex.glsl'),
                    fragmentShader: require('./shaders/Volume.lightmap.fragment.glsl'),
                    clipping: true,
                    depthWrite: false,
                    depthTest: false,
                }),
            );

            // for clipping planes
            quadRTT.applyMatrix4(modelMatrix);
            quadRTT.updateMatrixWorld();

            cameraRTT = new THREE.OrthographicCamera(
                lightMapRenderTargetSize / -2, lightMapRenderTargetSize / 2,
                lightMapRenderTargetSize / 2, lightMapRenderTargetSize / -2,
                -10000, 10000,
            );

            cameraRTT.position.z = 100;
            sceneRTT.add(quadRTT);

            object.refreshLightMap = function (direction) {
                if (textureRTT) {
                    const { renderer } = K3D.getWorld();
                    const cameraPosition = new THREE.Vector3();

                    if (direction) {
                        quadRTT.material.uniforms.lightDirection.value.fromArray(direction).normalize();
                    } else {
                        K3D.getWorld().camera.getWorldPosition(cameraPosition);
                        quadRTT.material.uniforms.lightDirection.value.copy(
                            translation.clone().sub(cameraPosition).normalize(),
                        );
                    }

                    K3D.getWorld().camera.updateMatrixWorld();

                    quadRTT.material.uniformsNeedUpdate = true;

                    renderer.clippingPlanes = [];
                    K3D.parameters.clippingPlanes.forEach((plane) => {
                        renderer.clippingPlanes.push(new THREE.Plane(new THREE.Vector3().fromArray(plane), plane[3]));
                    });

                    renderer.setRenderTarget(textureRTT);
                    renderer.clear(true, true, true);
                    renderer.render(sceneRTT, cameraRTT);
                    renderer.setRenderTarget(null);
                }
            };

            if (config.shadow === 'dynamic') {
                listenersId = K3D.on(K3D.events.BEFORE_RENDER, () => {
                    const now = new Date().getTime();

                    if (timeoutId) {
                        clearTimeout(timeoutId);
                    }

                    // check if we should updated shadow map because user interaction
                    if (now - lastShadowMapUpdated >= config.shadow_delay) {
                        object.refreshLightMap();
                        lastShadowMapUpdated = now;
                    } else {
                        // handle last update on end of user interaction
                        timeoutId = setTimeout(() => {
                            object.refreshLightMap();
                            lastShadowMapUpdated = now;
                            K3D.render();
                        }, Math.max(config.shadow_delay, 500));
                    }
                });
            }

            object.quadRTT = quadRTT;
            object.refreshLightMap();
        }

        object.onRemove = function () {
            if (quadRTT) {
                quadRTT.material.uniforms.volumeTexture.value.dispose();
                quadRTT.material.uniforms.volumeTexture.value = undefined;
            }

            object.material.uniforms.volumeTexture.value = undefined;
            object.material.uniforms.colormap.value.dispose();
            object.material.uniforms.colormap.value = undefined;
            jitterTexture.dispose();
            jitterTexture = undefined;

            if (sceneRTT) {
                sceneRTT = undefined;
            }

            if (cameraRTT) {
                cameraRTT = undefined;
            }

            if (textureRTT) {
                textureRTT.dispose();
                textureRTT = undefined;
            }

            if (listenersId) {
                K3D.off(K3D.events.BEFORE_RENDER, listenersId);
            }
        };

        return Promise.resolve(object);
    },

    update(config, changes, obj) {
        const resolvedChanges = {};

        if (typeof (changes.color_range) !== 'undefined' && !changes.color_range.timeSeries) {
            obj.material.uniforms.low.value = changes.color_range[0];
            obj.material.uniforms.high.value = changes.color_range[1];

            resolvedChanges.color_range = null;
        }

        if (typeof (changes.focal_length) !== 'undefined' && !changes.focal_length.timeSeries) {
            if ((obj.material.uniforms.focal_length.value === 0.0 && changes.focal_length !== 0.0)
                || changes.focal_length === 0.0) {
                // shader needs to be recompile
                return false;
            }
        }

        if (typeof (changes.volume) !== 'undefined' && !changes.volume.timeSeries) {
            if (obj.material.uniforms.volumeTexture.value.image.data.constructor === changes.volume.data.constructor
                && obj.material.uniforms.volumeTexture.value.image.width === changes.volume.shape[2]
                && obj.material.uniforms.volumeTexture.value.image.height === changes.volume.shape[1]
                && obj.material.uniforms.volumeTexture.value.image.depth === changes.volume.shape[0]) {
                obj.material.uniforms.volumeTexture.value.image.data = changes.volume.data;
                obj.material.uniforms.volumeTexture.value.needsUpdate = true;

                resolvedChanges.volume = null;
            }
        }

        if ((typeof (changes.color_map) !== 'undefined' && !changes.color_map.timeSeries)
            || (typeof (changes.opacity_function) !== 'undefined' && !changes.opacity_function.timeSeries)) {
            const canvas = colorMapHelper.createCanvasGradient(
                (changes.color_map && changes.color_map.data) || config.color_map.data,
                1024,
                (changes.opacity_function && changes.opacity_function.data) || config.opacity_function.data,
            );

            obj.material.uniforms.colormap.value.image = canvas;
            obj.material.uniforms.colormap.value.needsUpdate = true;

            if (obj.quadRTT) {
                obj.quadRTT.material.uniforms.colormap.value.image = canvas;
                obj.quadRTT.material.uniforms.colormap.value.needsUpdate = true;
            }

            resolvedChanges.color_map = null;
            resolvedChanges.opacity_function = null;
        }

        ['samples', 'alpha_coef', 'gradient_step', 'focal_plane', 'focal_length'].forEach((key) => {
            if (typeof (changes[key]) !== 'undefined' && !changes[key].timeSeries) {
                obj.material.uniforms[key].value = changes[key];
                resolvedChanges[key] = null;
            }
        });

        commonUpdate(config, changes, resolvedChanges, obj);

        if (areAllChangesResolve(changes, resolvedChanges)) {
            return Promise.resolve({ json: config, obj });
        }
        return false;
    },
};
